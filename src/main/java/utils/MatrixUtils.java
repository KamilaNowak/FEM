package utils;

import schemas.Element;
import schemas.Node;

import java.util.List;

public class MatrixUtils {
    private static final double pc = 0.57735026919;

    public static double[][] ksiEtaNetwork =
            {{-pc, -1, pc, -1},
                    {1, -pc, 1, pc},
                    {pc, 1, -pc, 1},
                    {-1, pc, -1, -pc}};


    public double calculateXJacobian(double[] derivatives, Element element) {
        double jacobianDerivativeX = 0;
        List<Node> elementNodes = element.getNodes();

        double[] x = new double[4];
        x[0] = elementNodes.get(0).getX();
        x[1] = elementNodes.get(1).getX();
        x[2] = elementNodes.get(2).getX();
        x[3] = elementNodes.get(3).getX();
        for (int i = 0; i < 4; i++)
            jacobianDerivativeX += derivatives[i] * x[i];
        return jacobianDerivativeX;
    }

    public double calculateYJacobian(double[] derivatives, Element element) {
        double jacobianDerivativeY = 0;
        List<Node> elementNodes = element.getNodes();

        double[] y = new double[4];
        y[0] = elementNodes.get(0).getY();
        y[1] = elementNodes.get(1).getY();
        y[2] = elementNodes.get(2).getY();
        y[3] = elementNodes.get(3).getY();
        for (int i = 0; i < 4; i++) jacobianDerivativeY += derivatives[i] * y[i];
        return jacobianDerivativeY;
    }
    public static double dNdx(double dYdKsi, double dYdEta, double dNdKsi, double dNdEta) {
        return dYdEta * dNdKsi + dYdKsi * dNdEta;
    }

    public static double dNdy(double dXdKsi, double dXdEta, double dNdKsi, double dNdEta) {
        return dXdEta * dNdKsi + dXdKsi * dNdEta;
    }

    public static double shapeFunctionPC1(double ksi, double ni) {
        return 0.25 * (1 - ksi) * (1 - ni);
    }

    public static double shapeFunctionPC2(double ksi, double ni) {
        return 0.25 * (1 + ksi) * (1 - ni);
    }

    public static double shapeFunctionPC3(double ksi, double ni) {
        return 0.25 * (1 + ksi) * (1 + ni);
    }

    public static double shapeFunctionPC4(double ksi, double ni) {
        return 0.25 * (1 - ksi) * (1 + ni);
    }


    public static double ksiDerivativePC1(double ni) {
        return -(1 - ni) / 4;
    }

    public static double ksiDerivativePC2(double ni) {
        return (1 - ni) / 4;
    }

    public static double ksiDerivativePC3(double ni) {
        return (1 + ni) / 4;
    }

    public static double ksiDerivativePC4(double ni) {
        return -(1 + ni) / 4;
    }

    public static double etaDerivativePC1(double ksi) {
        return -(1 - ksi) / 4;
    }

    public static double etaDerivativePC2(double ksi) {
        return -(1 + ksi) / 4;
    }

    public static double etaDerivativePC3(double ksi) {
        return (1 + ksi) / 4;
    }

    public static double etaDerivativePC4(double ksi) {
        return (1 - ksi) / 4;
    }
    public static double[] MatrixXVector(double M[][], double V[]){
        double [] newV = new double[V.length];
        for(int i = 0; i < M.length; i++){
            int value = 0;
            for(int j = 0; j < M[i].length; j++){
                value += M[i][j] * V[j];
            }
            newV[i] = value;
        }
        return newV;
    }
    public static double[][] multiplyMatrix2x4(double[][] A, double[][] B) {
        double[][] result = new double[2][4];

        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 4; j++)
                result[i][j] = 0;
        }

        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 4; j++) {
                for (int k = 0; k < 2; k++)
                    result[i][j] += (A[i][k] * B[k][j]);
            }
        }
        return result;
    }
    public static double[][] multiplyMatrix2x9(double[][] A, double[][] B) {
        double[][] result = new double[2][9];

        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 9; j++)
                result[i][j] = 0;
        }

        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 9; j++) {
                for (int k = 0; k < 2; k++)
                    result[i][j] += (A[i][k] * B[k][j]);
            }
        }
        return result;
    }

    public static double[][] multiplyMatrix4x4Transposition(double[] A, double[] B) {
        double[][] C = new double[4][4];

        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                C[i][j] = 0;
            }
        }

        for (int l = 0; l < 4; l++) {
            for (int i = 0; i < 4; i++) {
                C[l][i] += (A[l] * B[i]);
            }
        }

        return C;
    }
}

