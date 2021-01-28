package calculations;

import java.io.FileNotFoundException;

import readers.DataReader;
import org.la4j.LinearAlgebra;
import org.la4j.Matrix;
import org.la4j.Vector;
import org.la4j.inversion.MatrixInverter;
import schemas.Element;
import schemas.Network;
import schemas.Node;
import schemas.Temperatures;
import utils.matrixes.MatrixUtils;

import java.util.ArrayList;
import java.util.List;

import static java.lang.Math.pow;
import static java.lang.StrictMath.sqrt;

public class MatrixCalculates {

    public double[][] calculateJacobianForSingleElement(Element element, double[] ksi, double[] eta) {
        MatrixUtils matrixUtils = new MatrixUtils();
        double dXdKsi = matrixUtils.calculateXJacobian(ksi, element);
        double dXdEta = matrixUtils.calculateXJacobian(eta, element);
        double dYdKsi = matrixUtils.calculateYJacobian(ksi, element);
        double dYdEta = matrixUtils.calculateYJacobian(eta, element);

        // jakobian to macierz [2][2] o wyrazach dXDKsi... dYdEta.
        // Jakobian liczymy dla każdego węzła, czyli 4 dla jednego elementu, bo dla każdego punktu całkowania
        double[][] jacobian = new double[2][2];
        jacobian[0][0] = dXdKsi;
        jacobian[0][1] = dYdKsi;
        jacobian[1][0] = dXdEta;
        jacobian[1][1] = dYdEta;

        return jacobian;
    }

    // MACIERZ H
    public double[][] calculateHMatrix(Element element, double[][] dNdKsi, double[][] dNdEta, DataReader data) {
        double[][] H = new double[4][4];

        for (int index = 0; index < 4; index++) {

            double[][] dNdKsi_dNdEta = new double[2][4];

            //tworze jakobian dla pojedynczego elementu, [2][2]
            double[][] jacobi = calculateJacobianForSingleElement(element, dNdKsi[index], dNdEta[index]);
            double detJ = (jacobi[0][0] * jacobi[1][1]) - (jacobi[0][1] - jacobi[1][0]);

            //odwracam ten jakobian tak jak we wzorze i mnożę jego składniki przez 1/wyznacznik jakobianu
            double[][] reversedJacobi = new double[2][2];
            reversedJacobi[0][0] = jacobi[1][1];
            reversedJacobi[1][1] = jacobi[0][0];
            reversedJacobi[0][1] = -jacobi[0][1];
            reversedJacobi[1][0] = -jacobi[1][0];

            for (int l = 0; l < 2; l++) {
                for (int k = 0; k < 2; k++)
                    reversedJacobi[l][k] = reversedJacobi[l][k] * (1 / detJ);
            }
            //tu tworze macierz drugiego skladnika mnozenia dN_i/dKsi oraz dN_i/dEta
            for (int l = 0; l < 4; l++) {
                dNdKsi_dNdEta[0][l] = dNdKsi[index][l];
                dNdKsi_dNdEta[1][l] = dNdEta[index][l];
            }
            //mnoze dwa skladniki otrzymujac macierz dN_i/dx oraz dN_i/dy
            double[][] dNdXY = MatrixUtils.multiplyMatrix2x4(reversedJacobi, dNdKsi_dNdEta);
            double[] dNdX = new double[4];
            double[] dNdY = new double[4];

            //rodzielam wynik dN_i/dx oraz dN_i/dy na odzielne tablice
            for (int l = 0; l < 4; l++) {
                dNdX[l] = dNdXY[0][l];
                dNdY[l] = dNdXY[1][l];
            }
            //koncowe mnożenie macierzy wraz z transpozycją
            double[][] dNdX_dNdX_T = MatrixUtils.multiplyMatrix4x4Transposition(dNdX, dNdX);
            double[][] dNdY_dNdY_T = MatrixUtils.multiplyMatrix4x4Transposition(dNdY, dNdY);

            //macierz H to SUMA pojedynczych składników pomnożona przez współczynnik k i detJ
            for (int k = 0; k < 4; k++) {
                for (int j = 0; j < 4; j++)
                    H[k][j] += data.getK() * (dNdX_dNdX_T[k][j] + dNdY_dNdY_T[k][j]) * detJ;
                //k - wsp. przewodnictwa
            }
        }
        element.setH(H);
        return H;
    }

    // MACIERZ GLOBAL  H
    public double[][] calculateGlobalHMatrix(Network network) {

        int sizeOfIds = network.getElements().get(0).getIDs().size();
        int matrixSize = network.getNodes().size();
        double[][] H_global = new double[matrixSize][matrixSize];

        for (Element element : network.getElements()) {
            for (int i = 0; i < sizeOfIds; i++) {
                for (int j = 0; j < sizeOfIds; j++) {
                    H_global[element.getIDs().get(i) - 1][element.getIDs().get(j) - 1]
                            += (element.getH()[i][j] + element.getHbc()[i][j]);
                }
            }
        }
        return H_global;
    }

    // MACIERZ GLOBAL C
    public double[][] generateC(Element element, double[][] N, double[][] dNdKsi, double[][] dNdEta, DataReader data) {
        double[][] C = new double[4][4];

        for (int index = 0; index < 4; index++) {
            // C = c * {N} * {N_transp} * detJ
            double[][] jacobi = calculateJacobianForSingleElement(element, dNdKsi[index], dNdEta[index]);
            double detJ = (jacobi[0][0] * jacobi[1][1]) - (jacobi[0][1] - jacobi[1][0]);
            double[][] final_N = MatrixUtils.multiplyMatrix4x4Transposition(N[index], N[index]);

            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                    C[i][j] += final_N[i][j] * data.getC() * data.getRo() * detJ;
                }
            }
        }
        element.setC(C);
        return C;
    }

    public double[][] generateCFor3Points(Element element, double[][] N, double[][] dNdKsi, double[][] dNdEta, DataReader data) {
        double[][] C = new double[4][4];

        for (int index = 0; index < 9; index++) {

            double[][] jacobi = calculateJacobianForSingleElement(element, dNdKsi[index], dNdEta[index]);
            double detJ = (jacobi[0][0] * jacobi[1][1]) - (jacobi[0][1] - jacobi[1][0]);
            double[][] final_N = MatrixUtils.multiplyMatrix9x4Transposition(N[index], N[index]);

            for (int k = 0; k < 4; k++) {
                for (int j = 0; j < 4; j++)
                    if (index == 0 || index == 2 || index == 6 || index == 8) {
                        C[k][j] += (final_N[k][j] * data.getC() * data.getRo() * detJ)
                                * data.getWeights3().get(0) * data.getWeights3().get(2); // 0,55 - 0,55
                    } else if (index == 4) {
                        C[k][j] += (final_N[k][j] * data.getC() * data.getRo() * detJ)
                                * data.getWeights3().get(1) * data.getWeights3().get(1); // 0,88 - 0,88
                    } else {
                        C[k][j] += (final_N[k][j] * data.getC() * data.getRo() * detJ)
                                * data.getWeights3().get(0) * data.getWeights3().get(1);  // 0,55 - 0,88
                    }
            }
        }
        element.setC(C);
        return C;
    }

    // MATRIX GLOBAL C
    public double[][] calculateGlobalCMatrix(Network network) {
        int sizeOFIds = network.getElements().get(0).getIDs().size();
        int matrixSize = network.getNodes().size();

        double[][] C_global = new double[matrixSize][matrixSize];
        for (Element element : network.getElements()) {
            for (int i = 0; i < sizeOFIds; i++) {
                for (int j = 0; j < sizeOFIds; j++) {
                    C_global[element.getIDs().get(i) - 1][element.getIDs().get(j) - 1] += element.getC()[i][j];
                }
            }
        }
        return C_global;
    }

    // MATRIX  HBC
    public double[][] calculateHbcMatrix(Element element, DataReader data) throws FileNotFoundException {
        List<Node> nodeList = element.getNodes();
        double[][] Hbc = new double[4][4];

        boolean[] bcSurface = bc(nodeList);

        // Hbc = alfa * {N} * {N_transp} * detS
        for (int i = 0; i < 4; i++) {
            if (bcSurface[i] == true) {
                double[] ksiEtaSurface = MatrixUtils.ksiEtaSurface2p()[i]; //biore pkt calkowania ze ścianki
                double[] ksi = {ksiEtaSurface[0], ksiEtaSurface[2]};
                double[] eta = {ksiEtaSurface[1], ksiEtaSurface[3]};

                //macierz N [2][4] licze funkcje kształtu
                double[][] surface = Interpolation.calculateShapeFunctionsSurface(ksi, eta, data);

                //obliczanie l (tw. pitagorasa)
                double sideLength = getSideLength(nodeList, i);
                double detS = sideLength / 2;

                double[][] side0 = MatrixUtils.multiplyMatrix4x4Transposition(surface[0], surface[0]);
                double[][] side1 = MatrixUtils.multiplyMatrix4x4Transposition(surface[1], surface[1]);

                // sumowanie macierzy
                double[][] sidesSum = MatrixUtils.sumMatrix4x4(side0, side1);

                for (int l = 0; l < sidesSum.length; l++) {
                    for (int m = 0; m < sidesSum[i].length; m++)
                        Hbc[l][m] += (sidesSum[l][m] * detS * data.getAlfa());
                }
            }
        }
        element.setHbc(Hbc);
        return Hbc;
    }

    // MATRIX P
    public double[] calculatePMatrix(Element element, DataReader data) throws FileNotFoundException {
        List<Node> nodeList = element.getNodes();
        double[][] P = new double[2][(int) data.getNW()];
        double[] Pvector = new double[(int) data.getNW()];

        boolean[] bcSurface = bc(nodeList);
        for (int i = 0; i < 4; i++) {
            if (bcSurface[i] == true) {

                double[] ksiEtaSurface = MatrixUtils.ksiEtaSurface2p()[i];
                double[] ksi = {ksiEtaSurface[0], ksiEtaSurface[2]};
                double[] eta = {ksiEtaSurface[1], ksiEtaSurface[3]};

                double[][] surface = Interpolation.calculateShapeFunctionsSurface(ksi, eta, data);

                double sideLength = getSideLength(nodeList, i);
                double detS = sideLength / 2;

                for (int l = 0; l < surface.length; l++) {
                    for (int m = 0; m < surface[l].length; m++)
                        P[l][m] += (-1) * data.getAlfaForP() * surface[l][m] * data.getTempAlfa() * detS;

                }
            }
        }
        //najpierw mnoze wszystko na raz, potem wyłuskuje sume wierszy do wektora P
        for (int index = 0; index < 4; index++)
            Pvector[index] = (P[0][index] + P[1][index]); // sumuję jak we wzorze w1 + w2

        element.setP(Pvector);

        return Pvector;
    }

    public Object[] aggregation(DataReader data, Network network, double[][] H_global, double[][] C_global) {

        int SIZE = (int) (data.getNH() * data.getNW());
        double[] P_Global = new double[SIZE];
        double[][] H = new double[SIZE][SIZE];
        double[][] C = new double[SIZE][SIZE];
        double[] temp = new double[SIZE];

        // liczę P Globalne
        for (Element element : network.getElements()) {
            List<Integer> ID = element.getIDs();
            for (int k = 0; k < 4; k++)
                P_Global[ID.get(k) - 1] += element.getP()[k];
        }
       //pobieram temperastury do tablicy
        for (int i = 0; i < network.getNodes().size(); i++)
            temp[i] = network.getNodes().get(i).getT();

        // ([H] + [C]/dt) -([C]/dt) * {t} +{P}=0
        for (int i = 0; i < SIZE; i++) {
            for (int j = 0; j < SIZE; j++) {
                H[i][j] = H_global[i][j] + (C_global[i][j] / data.getSimulationStepTime()); // składnik ([H] + [C]/dt)
                C[i][j] = (C_global[i][j] / data.getSimulationStepTime());  // składnik [C]/dt
            }
        }

        double[] CxT0 = MatrixUtils.MatrixXVector(C, temp); // składnik ([C]/dt) * {t}
        double[] P_Final = new double[SIZE];

        for (int i = 0; i < SIZE; i++)
            P_Final[i] = (-1) * CxT0[i] + P_Global[i]; // {P^} = {P} - ([C]/ dTao) * {t}

        System.out.println("Max and min temperature in each step");
        List<Temperatures> temperatures = calculateTemperatures(H, C, P_Global, network, temp, data);

        for (Temperatures t : temperatures)
            System.out.println(t.getTime() + "\t\t" + t.getMinTemperature() + "\t\t" + t.getMaxTemperature());
        System.out.println();

        return new Object[]{H, P_Final, C};
    }

    public List<Temperatures> calculateTemperatures(double[][] H, double[][] C, double[] P_Global, Network network, double[] temp, DataReader data) {
        System.out.println("Time[s]\t\t\tMinTemp[s]\t\t\t\tMaxTemp[s]");

        List<Temperatures> temperatures = new ArrayList<>();

        for (int i = 0; i < data.getSimulationTime() / data.getSimulationStepTime(); i++) {

            for (int j = 0; j < network.getNodes().size(); j++) temp[j] = network.getNodes().get(j).getT();
            //([H + [C]/dtao) *{t1} - ([C]/dTao)* {t0} + {P}
            Vector temperature = Vector.fromArray(temp);
            Vector P = Vector.fromArray(P_Global);
            Matrix H_tmp = Matrix.from2DArray(H);
            Matrix C_tmp = Matrix.from2DArray(C);

            Vector Cxt0 = C_tmp.multiply(temperature);
            Cxt0 = Cxt0.subtract(P);

            MatrixInverter inverter = H_tmp.withInverter(LinearAlgebra.GAUSS_JORDAN);
            Matrix H_inversed = inverter.inverse();
            Vector result = Cxt0.multiply(H_inversed); // wynik - wektor temperatur {t1}

            for (int j = 0; j < network.getNodes().size(); j++) network.getNodes().get(j).setT(result.get(j));

            double time = data.getSimulationStepTime() * (i + 1);

            double minTemp = result.min();
            double maxTemp = result.max();
            temperatures.add(new Temperatures(time, minTemp, maxTemp));
        }
        return temperatures;
    }

    double getSideLength(List<Node> nodeList, int index) {
        if (index == 3) { // zeby nie wyjsc poza zakres tablicy bo rzuci wyjątek, więc biore ostatni i pierwszy
            return sqrt(
                    pow(nodeList.get(0).getX() - nodeList.get(3).getX(), 2)
                            + pow(nodeList.get(0).getY() - nodeList.get(3).getY(), 2)
            );
        } else return sqrt(
                pow(nodeList.get(index + 1).getX() - nodeList.get(index).getX(), 2)
                        + pow(nodeList.get(index + 1).getY() - nodeList.get(index).getY(), 2)
        );
    }

    public boolean[] bc(List<Node> nodeList) {
        boolean[] bcSurface = new boolean[4];
        boolean[] bc = new boolean[4];

        for (int i = 0; i < nodeList.size(); i++) {
            bc[i] = nodeList.get(i).isBc(); // sprawdzam czy jest to node na sciance
            System.out.print(" - " + bc[i]);
        }

        for (int i = 0; i <= 3; i++) {
            if (i == 3) {
                if (bc[0] == true && bc[3] == true) {
                    bcSurface[3] = true;
                } // czy scianka
                // zeby nie wyjsc poza zakres tablicy bo rzuci wyjątek, więc biore ostatni i pierwszy
            } else if (bc[i] == true && bc[i + 1] == true) {
                bcSurface[i] = true;
            }
        }
        return bcSurface;
    }

    public double[][] calculateHMatrix3Points(Element element, double[][] dNdKsi, double[][] dNdEta, DataReader data) {
        double[][] H = new double[4][4];

        for (int index = 0; index < 9; index++) {

            double[][] dNdKsi_dNdEta = new double[2][4];
            double[] dNdX = new double[4];
            double[] dNdY = new double[4];

            double[][] jacobi = calculateJacobianForSingleElement(element, dNdKsi[index], dNdEta[index]);

            double detJ = (jacobi[0][0] * jacobi[1][1]) - (jacobi[0][1] - jacobi[1][0]);
            double[][] reversedJacobi = new double[2][2];

            reversedJacobi[0][0] = jacobi[1][1];
            reversedJacobi[1][1] = jacobi[0][0];

            reversedJacobi[0][1] = -jacobi[0][1];
            reversedJacobi[1][0] = -jacobi[1][0];

            for (int l = 0; l < 2; l++) {
                for (int k = 0; k < 2; k++)
                    reversedJacobi[l][k] = reversedJacobi[l][k] * (1 / detJ);
            }

            for (int l = 0; l < 4; l++) {
                dNdKsi_dNdEta[0][l] = dNdKsi[index][l];
                dNdKsi_dNdEta[1][l] = dNdEta[index][l];
            }
            double[][] dNdXY = MatrixUtils.multiplyMatrix2x4(reversedJacobi, dNdKsi_dNdEta);

            for (int l = 0; l < 4; l++) {
                dNdX[l] = dNdXY[0][l];
                dNdY[l] = dNdXY[1][l];
            }

            double[][] X = MatrixUtils.multiplyMatrix4x4Transposition(dNdX, dNdX);
            double[][] Y = MatrixUtils.multiplyMatrix4x4Transposition(dNdY, dNdY);

            // Mnożenie z wagami 0,555, 0,888
            for (int k = 0; k < 4; k++) {
                for (int j = 0; j < 4; j++)
                    if (index == 0 || index == 2 || index == 6 || index == 8) {
                        H[k][j] += (data.getK() * (X[k][j] + Y[k][j]) * detJ)
                                * data.getWeights3().get(0) * data.getWeights3().get(2); // 0,55 - 0,55
                    } else if (index == 4) {
                        H[k][j] += (data.getK() * (X[k][j] + Y[k][j]) * detJ)
                                * data.getWeights3().get(1) * data.getWeights3().get(1); // 0,88 - 0,88
                    } else {
                        H[k][j] += (data.getK() * (X[k][j] + Y[k][j]) * detJ)
                                * data.getWeights3().get(0) * data.getWeights3().get(1);  // 0,55 - 0,88
                    }
            }
        }
        element.setH(H);
        return H;
    }
}