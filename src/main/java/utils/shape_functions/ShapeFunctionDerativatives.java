package utils.shape_functions;

public class ShapeFunctionDerativatives {

    public static final double INIT_VALUE = 0.25;

    //pochodne funkcji kształtu ksi i eta dla kazdego punktu calkowania
    public static double ksiDerivativePC1(double eta) { return -(1 - eta) * INIT_VALUE; }

    public static double etaDerivativePC1(double ksi) { return -(1 - ksi) * INIT_VALUE; }

    public static double ksiDerivativePC2(double eta) {
        return (1 - eta) * INIT_VALUE;
    }

    public static double ksiDerivativePC3(double eta) {
        return (1 + eta) * INIT_VALUE;
    }

    public static double ksiDerivativePC4(double eta) {
        return -(1 + eta) * INIT_VALUE;
    }



    public static double etaDerivativePC2(double ksi) {
        return -(1 + ksi) * INIT_VALUE;
    }

    public static double etaDerivativePC3(double ksi) {
        return (1 + ksi) * INIT_VALUE;
    }

    public static double etaDerivativePC4(double ksi) {
        return (1 - ksi) * INIT_VALUE;
    }
}
