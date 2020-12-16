package utils;

public class ShapeFunctions {

    public static double PC1(double ksi, double eta) { return 0.25 * (1 - ksi) * (1 - eta); }

    public static double PC2(double ksi, double eta) {
        return 0.25 * (1 + ksi) * (1 - eta);
    }

    public static double PC3(double ksi, double eta) {
        return 0.25 * (1 + ksi) * (1 + eta);
    }

    public static double PC4(double ksi, double eta) {
        return 0.25 * (1 - ksi) * (1 + eta);
    }
}
