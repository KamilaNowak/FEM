package printable;

import java.text.DecimalFormat;

public class Print {

    public static void printMatrixNxM(double[][] m, int N, int M) {
        DecimalFormat df = new DecimalFormat("#.##");
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < M; j++)
                System.out.printf("%9.8s" ,df.format(m[i][j]) + "\t");
            System.out.println();
        }
    }

    public static void printArray(double[] V, int l) {
        DecimalFormat df = new DecimalFormat("#.##");
        for (int i = 0; i < l; i++)
            System.out.printf("%9.9s", df.format(V[i]) + "\t");
        System.out.println();
    }
}
