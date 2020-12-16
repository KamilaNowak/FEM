package utils;

public class Derativatives {

    public static double ksiDerivativePC1(double eta) {
        return -(1 - eta) / 4;
    }

    public static double ksiDerivativePC2(double eta) {
        return (1 - eta) / 4;
    }

    public static double ksiDerivativePC3(double eta) {
        return (1 + eta) / 4;
    }

    public static double ksiDerivativePC4(double eta) {
        return -(1 + eta) / 4;
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
}
