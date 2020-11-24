package calculations;

import input.DataReader;
import printable.Print;
import schemas.Element;
import schemas.Network;
import sun.text.resources.cldr.pa.FormatData_pa_Arab;
import utils.MatrixUtils;

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
        int sizeOdIds = network.getElements().get(0).getiDs().size();
        int matrixSize = network.getNodes().size();

        double[][] H_global = new double[matrixSize][matrixSize];
        for (Element element : network.getElements()) {
            for (int i = 0; i < sizeOdIds; i++) {
                for (int j = 0; j < sizeOdIds; j++) {
                    H_global[element.getiDs().get(i) - 1][element.getiDs().get(j) - 1] += element.getH()[i][j];
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
            double[][] final_N = MatrixUtils.multiplyMatrix4x4Transposition(N[index], N[index]);

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
        int sizeOdIds = network.getElements().get(0).getiDs().size();
        int matrixSize = network.getNodes().size();

        double[][] C_global = new double[matrixSize][matrixSize];
        for (Element element : network.getElements()) {
            for (int i = 0; i < sizeOdIds; i++) {
                for (int j = 0; j < sizeOdIds; j++) {
                    C_global[element.getiDs().get(i) - 1][element.getiDs().get(j) - 1] += element.getC()[i][j];
                }
            }
        }
        return C_global;
    }


}