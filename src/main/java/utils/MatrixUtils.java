package utils;

import input.DataReader;
import schemas.Element;
import schemas.Node;
import java.io.FileNotFoundException;
import java.util.List;

public class MatrixUtils {

    private static DataReader data;

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


    public static double[][] ksiEtaSurface2p() throws FileNotFoundException {
        data = new DataReader();
        double pc1 = data.getPc2().get(0); //-0.57735026919;
        double pc2 = data.getPc2().get(1); //0.57735026919;

        return new double[][]{ //
                {pc1, -1, pc2, -1},
                {1, pc1, 1, pc2},
                {pc2, 1, pc1, 1},
                {-1, pc2, -1, pc1}
        };
    }

    public static double[][] sumMatrix4x4(double[][] arr1, double[][] arr2) {
        double[][] sumArray = new double[4][4];
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++)
                sumArray[i][j] = arr1[i][j] + arr2[i][j];
        }
        return sumArray;
    }

    public static double[] MatrixXVector(double[][] matrix, double[] vector) {
        double[] resultVector = new double[vector.length];
        for (int i = 0; i < matrix.length; i++) {
            int value = 0;
            for (int j = 0; j < matrix[i].length; j++)
                value += matrix[i][j] * vector[j];
            resultVector[i] = value;
        }
        return resultVector;
    }

    public static double[][] multiplyMatrix2x4(double[][] arr1, double[][] arr2) {
        double[][] resultMatrix = new double[2][4];
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 4; j++)
                resultMatrix[i][j] = 0;
        }
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 4; j++) {
                for (int k = 0; k < 2; k++)
                    resultMatrix[i][j] += (arr1[i][k] * arr2[k][j]);
            }
        }
        return resultMatrix;
    }

    public static double[][] multiplyMatrix2x9(double[][] arr1, double[][] arr2) {
        double[][] resultMatrix = new double[2][9];
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 9; j++)
                resultMatrix[i][j] = 0;
        }
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 9; j++) {
                for (int k = 0; k < 2; k++)
                    resultMatrix[i][j] += (arr1[i][k] * arr2[k][j]);
            }
        }
        return resultMatrix;
    }

    public static double[][] multiplyMatrix4x4Transposition(double[] arr1, double[] arr2) {
        double[][] resultArray = new double[4][4];
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++)
                resultArray[i][j] = 0;
        }
        for (int l = 0; l < 4; l++) {
            for (int i = 0; i < 4; i++)
                resultArray[l][i] += (arr1[l] * arr2[i]);
        }
        return resultArray;
    }

    public static double[][] multiplyMatrix9x4Transposition(double[] arr1, double[] arr2) {
        double[][] resultArray = new double[9][4];
        for (int i = 0; i < 9; i++) {
            for (int j = 0; j < 4; j++)
                resultArray[i][j] = 0;
        }
        for (int l = 0; l < 9; l++) {
            for (int i = 0; i < 4; i++)
                resultArray[l][i] += (arr1[i] * arr2[i]);
        }
        return resultArray;
    }

}

