package printable;

import java.text.DecimalFormat;

public class Print {

    public static void printMatrixNxM(double[][] matrix, int len1, int len2) {
        DecimalFormat df = new DecimalFormat("#.##");
        for (int i = 0; i < len1; i++) {
            for (int j = 0; j < len2; j++)
                System.out.printf("%9.8s" ,df.format(matrix[i][j]) + "\t");
            System.out.println();
        }
    }

    public static void printVec(double[] vector) {
        DecimalFormat df = new DecimalFormat("#.##");
        for (double v : vector) System.out.printf("%8.9s", df.format(v) + "\t");
        System.out.println();
    }

    public static void printVec(double[] vector, int length) {
        for (int i = 0; i < length; i++) {
            System.out.print(vector[i] +"\t");
        }
        System.out.println();
    }
}
