package calculations;

import java.io.FileNotFoundException;

import input.DataReader;
import javafx.util.Pair;
import org.la4j.LinearAlgebra;
import org.la4j.Matrices;
import org.la4j.Matrix;
import org.la4j.Vector;
import org.la4j.inversion.MatrixInverter;
import printable.Print;
import schemas.Element;
import schemas.Network;
import schemas.Node;
import schemas.Temperatures;
import utils.MatrixUtils;

import java.text.DecimalFormat;
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

        double[][] jacobian = new double[2][2];
        jacobian[0][0] = dXdKsi;
        jacobian[0][1] = dYdKsi;
        jacobian[1][0] = dXdEta;
        jacobian[1][1] = dYdEta;

        return jacobian;
    }

    public double[][] calculateJacobianForSingleElement3Points(Element element, double[] ksi, double[] eta) {

        MatrixUtils matrixUtils = new MatrixUtils();
        double[][] jacobi = new double[2][2];

        double dXdKsi = matrixUtils.calculateXJacobian(ksi, element);
        double dXdEta = matrixUtils.calculateXJacobian(eta, element);
        double dYdKsi = matrixUtils.calculateYJacobian(ksi, element);
        double dYdEta = matrixUtils.calculateYJacobian(eta, element);

        jacobi[0][0] = dXdKsi;
        jacobi[0][1] = dYdKsi;
        jacobi[1][0] = dXdEta;
        jacobi[1][1] = dYdEta;

        return jacobi;
    }

    public double[][] calculateHMatrix(Element element, double[][] dNdKsi, double[][] dNdEta, DataReader data) {
        double[][] H = new double[4][4];

        for (int index = 0; index < 4; index++) {

            double[][] dNdKsi_dNdEta = new double[2][4];
            double[][] jacobi = calculateJacobianForSingleElement(element, dNdKsi[index], dNdEta[index]);

            double detJ = (jacobi[0][0] * jacobi[1][1]) - (jacobi[0][1] - jacobi[1][0]);

            double[][] reversedJacobi = new double[2][2];
            reversedJacobi[0][0] = jacobi[1][1];
            reversedJacobi[0][1] = -jacobi[0][1];
            reversedJacobi[1][0] = -jacobi[1][0];
            reversedJacobi[1][1] = jacobi[0][0];

            reversedJacobi[0][0] = reversedJacobi[0][0] * (1 / detJ);
            reversedJacobi[0][1] = reversedJacobi[0][1] * (1 / detJ);
            reversedJacobi[1][0] = reversedJacobi[1][0] * (1 / detJ);
            reversedJacobi[1][1] = reversedJacobi[1][1] * (1 / detJ);

            dNdKsi_dNdEta[0][0] = dNdKsi[index][0];
            dNdKsi_dNdEta[0][1] = dNdKsi[index][1];
            dNdKsi_dNdEta[0][2] = dNdKsi[index][2];
            dNdKsi_dNdEta[0][3] = dNdKsi[index][3];

            dNdKsi_dNdEta[1][0] = dNdEta[index][0];
            dNdKsi_dNdEta[1][1] = dNdEta[index][1];
            dNdKsi_dNdEta[1][2] = dNdEta[index][2];
            dNdKsi_dNdEta[1][3] = dNdEta[index][3];

            double[][] dNdXY = MatrixUtils.multiplyMatrix2x4(reversedJacobi, dNdKsi_dNdEta);

            double[] dNdX = new double[4];
            double[] dNdY = new double[4];

            dNdX[0] = dNdXY[0][0];
            dNdX[1] = dNdXY[0][1];
            dNdX[2] = dNdXY[0][2];
            dNdX[3] = dNdXY[0][3];

            dNdY[0] = dNdXY[1][0];
            dNdY[1] = dNdXY[1][1];
            dNdY[2] = dNdXY[1][2];
            dNdY[3] = dNdXY[1][3];

            double[][] X = MatrixUtils.multiplyMatrix4x4Transposition(dNdX, dNdX);
            double[][] Y = MatrixUtils.multiplyMatrix4x4Transposition(dNdY, dNdY);

            for (int k = 0; k < 4; k++) {
                for (int j = 0; j < 4; j++)
                    H[k][j] += data.getK() * (X[k][j] + Y[k][j]) * detJ;
            }
        }
        element.setH(H);
        return H;
    }

    public double[][] calculateHMatrix3Points(Element element, double[][] dNdKsi, double[][] dNdEta, DataReader data) {
        double[][] H = new double[4][4];

        for (int index = 0; index < 9; index++) {

            double[][] dNdKsi_dNdEta = new double[2][4];
            double[] dNdX = new double[4];
            double[] dNdY = new double[4];

            double[][] jacobi = calculateJacobianForSingleElement3Points(element, dNdKsi[index], dNdEta[index]);

            double detJ = (jacobi[0][0] * jacobi[1][1]) - (jacobi[0][1] - jacobi[1][0]);

            double[][] reversedJacobi = new double[2][2];
            reversedJacobi[0][0] = jacobi[1][1];
            reversedJacobi[0][1] = -jacobi[0][1];
            reversedJacobi[1][0] = -jacobi[1][0];
            reversedJacobi[1][1] = jacobi[0][0];

            reversedJacobi[0][0] = reversedJacobi[0][0] * (1 / detJ);
            reversedJacobi[0][1] = reversedJacobi[0][1] * (1 / detJ);
            reversedJacobi[1][0] = reversedJacobi[1][0] * (1 / detJ);
            reversedJacobi[1][1] = reversedJacobi[1][1] * (1 / detJ);

            dNdKsi_dNdEta[0][0] = dNdKsi[index][0];
            dNdKsi_dNdEta[0][1] = dNdKsi[index][1];
            dNdKsi_dNdEta[0][2] = dNdKsi[index][2];
            dNdKsi_dNdEta[0][3] = dNdKsi[index][3];

            dNdKsi_dNdEta[1][0] = dNdEta[index][0];
            dNdKsi_dNdEta[1][1] = dNdEta[index][1];
            dNdKsi_dNdEta[1][2] = dNdEta[index][2];
            dNdKsi_dNdEta[1][3] = dNdEta[index][3];

            double[][] dNdXY = MatrixUtils.multiplyMatrix2x4(reversedJacobi, dNdKsi_dNdEta);

            dNdX[0] = dNdXY[0][0];
            dNdX[1] = dNdXY[0][1];
            dNdX[2] = dNdXY[0][2];
            dNdX[3] = dNdXY[0][3];

            dNdY[0] = dNdXY[1][0];
            dNdY[1] = dNdXY[1][1];
            dNdY[2] = dNdXY[1][2];
            dNdY[3] = dNdXY[1][3];

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

    public double[][] calculateGlobalHMatrix(Network network) {
        int sizeOdIds = network.getElements().get(0).getIDs().size();
        int matrixSize = network.getNodes().size();

        double[][] H_global = new double[matrixSize][matrixSize];
        for (Element element : network.getElements()) {
            for (int i = 0; i < sizeOdIds; i++) {
                for (int j = 0; j < sizeOdIds; j++) {
                    H_global[element.getIDs().get(i) - 1][element.getIDs().get(j) - 1] += (element.getH()[i][j] + element.getHbc()[i][j]);
                }
            }
        }
        return H_global;
    }

    ////// C Matrix  //////
    public double[][] generateC(Element element, double[][] N, double[][] dNdKsi, double[][] dNdEta, DataReader data) {
        double[][] C = new double[4][4];

        for (int index = 0; index < 4; index++) {

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

    public double[][] calculateGlobalCMatrix(Network network) {
        int sizeOdIds = network.getElements().get(0).getIDs().size();
        int matrixSize = network.getNodes().size();

        double[][] C_global = new double[matrixSize][matrixSize];
        for (Element element : network.getElements()) {
            for (int i = 0; i < sizeOdIds; i++) {
                for (int j = 0; j < sizeOdIds; j++) {
                    C_global[element.getIDs().get(i) - 1][element.getIDs().get(j) - 1] += element.getC()[i][j];
                }
            }
        }
        return C_global;
    }

    public boolean[] setBcSurface(List<Node> nodeList) {
        boolean[] bcSurface = new boolean[4];
        boolean[] bc = new boolean[4];

        for (int i = 0; i < nodeList.size(); i++)
            bc[i] = nodeList.get(i).isBc();

        for (int i = 0; i <= 3; i++) {
            if (i == 3) {
                if (bc[0] && bc[3])
                    bcSurface[3] = true;
            } else if (bc[i] && bc[i + 1]) {
                bcSurface[i] = true;
            }
        }
        return bcSurface;
    }

    public double[][] calculateHbcMatrix(Element element, DataReader data) throws FileNotFoundException {
        List<Node> nodeList = element.getNodes();
        double[][] Hbc = new double[4][4];

        boolean[] bcSurface = setBcSurface(nodeList);

        for (int i = 0; i < 4; i++) {
            if (bcSurface[i] == true) {
                double[] ksiEtaSurface = MatrixUtils.ksiEtaSurface2p()[i];
                double[] ksi = {ksiEtaSurface[0], ksiEtaSurface[2]};
                double[] eta = {ksiEtaSurface[1], ksiEtaSurface[3]};

                //macierz N [2][4]
                double[][] surface = Interpolation.calculateShapeFunctionsSurface(ksi, eta);

                //obliczanie l
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

    double getSideLength(List<Node> nodeList, int index) {
        if (index == 3) {
            return sqrt(
                    pow(nodeList.get(0).getX() - nodeList.get(3).getX(), 2)
                            + pow(nodeList.get(0).getY() - nodeList.get(3).getY(), 2)
            );
        } else return sqrt(
                pow(nodeList.get(index + 1).getX() - nodeList.get(index).getX(), 2)
                        + pow(nodeList.get(index + 1).getY() - nodeList.get(index).getY(), 2)
        );
    }


    public double[] calculatePMatrix(Element element, DataReader data) throws FileNotFoundException {
        List<Node> nodeList = element.getNodes();
        double[][] P = new double[2][(int) data.getNW()];
        double[] Pvector = new double[(int) data.getNW()];

        boolean[] bcSurface = setBcSurface(nodeList);
        for (int i = 0; i < 4; i++) {
            if (bcSurface[i] == true) {

                double[] ksiEtaSurface = MatrixUtils.ksiEtaSurface2p()[i];
                double[] ksi = {ksiEtaSurface[0], ksiEtaSurface[2]};
                double[] eta = {ksiEtaSurface[1], ksiEtaSurface[3]};

                double[][] surface = Interpolation.calculateShapeFunctionsSurface(ksi, eta);

                double sideLength = getSideLength(nodeList, i);
                double detS = sideLength / 2;

                for (int l = 0; l < surface.length; l++) {
                    for (int m = 0; m < surface[l].length; m++)
                        P[l][m] += (-1) * data.getAlfaForP() * surface[l][m] * data.getTempAlfa() * detS;
                }
            }
        }

        for (int index = 0; index < 4; index++) {
            Pvector[index] = (P[0][index] + P[1][index]);
        }
        element.setP(Pvector);

        return Pvector;
    }

    public Object[] aggregation(DataReader data, Network network, double[][] H_global, double[][] C_global) {

        int SIZE = (int) (data.getNH() * data.getNW());
        double[] P_Global = new double[SIZE];
        double[][] H = new double[SIZE][SIZE];
        double[][] C = new double[SIZE][SIZE];
        double[] temp = new double[SIZE];

        Pair<double[][],  double[]> result = new Pair(null, null);

        for (Element element : network.getElements()) {
            List<Integer> ID = element.getIDs();
            for (int k = 0; k < 4; k++)
                P_Global[ID.get(k) - 1] += element.getP()[k];
        }

        for (int i = 0; i < network.getNodes().size(); i++)
            temp[i] = network.getNodes().get(i).getT();

        for (int i = 0; i < SIZE; i++) {
            for (int j = 0; j < SIZE; j++) {
                H[i][j] = H_global[i][j] + (C_global[i][j] / data.getSimulationStepTime());
                C[i][j] = (C_global[i][j] / data.getSimulationStepTime());
            }
        }

        double[] CxT0 = MatrixUtils.MatrixXVector(C, temp);
        double[] P_Final = new double[SIZE];

        System.out.println("C * t0");
        for (int i = 0; i < SIZE; i++) {
            P_Final[i] = CxT0[i] - P_Global[i];
        }

        System.out.println("Max and min temperature in each step");
        List<Temperatures> temperatures = calculateTemperatures(H, C, P_Global, network, temp, data);
        for (Temperatures t : temperatures)
            System.out.println(t.getTime() + "\t\t" + t.getMinTemperature() + "\t\t" + t.getMaxTemperature());
        System.out.println();
        return new Object[]{H, P_Final,C};
    }

    public List<Temperatures> calculateTemperatures(double[][] H, double[][] C, double[] P_Global, Network network, double[] temp, DataReader data) {
        System.out.println("Time[s]\t\t\tMinTemp[s]\t\t\t\tMaxTemp[s]");

        List<Temperatures> temperatures = new ArrayList<>();

        for (int i = 0; i < data.getSimulationTime() / data.getSimulationStepTime(); i++) {

            for (int j = 0; j < network.getNodes().size(); j++) temp[j] = network.getNodes().get(j).getT();

            Vector temperature = Vector.fromArray(temp);
            Vector P = Vector.fromArray(P_Global);
            Matrix H_tmp = Matrix.from2DArray(H);
            Matrix C_tmp = Matrix.from2DArray(C);

            Vector Cxt0 = C_tmp.multiply(temperature);
            Cxt0 = Cxt0.subtract(P);

            MatrixInverter inverter = H_tmp.withInverter(LinearAlgebra.GAUSS_JORDAN);
            Matrix H_inversed = inverter.inverse();
            Vector result = Cxt0.multiply(H_inversed);

            for (int j = 0; j < network.getNodes().size(); j++) network.getNodes().get(j).setT(result.get(j));

            double time = data.getSimulationStepTime() * (i + 1);
            double minTemp = result.min();
            double maxTemp = result.max();
            temperatures.add(new Temperatures(time, minTemp, maxTemp));
        }
        return temperatures;
    }
}