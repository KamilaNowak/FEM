package utils.shape_functions;

public class ShapeFunctions {

    //funkcje kształtu dla kazdego punktu całkowaniw
    public static final double INIT_VALUE = 0.25;

    public static double PC1(double ksi, double eta) { return INIT_VALUE * (1 - ksi) * (1 - eta); }

    public static double PC2(double ksi, double eta) {
        return INIT_VALUE * (1 + ksi) * (1 - eta);
    }

    public static double PC3(double ksi, double eta) {
        return INIT_VALUE * (1 + ksi) * (1 + eta);
    }

    public static double PC4(double ksi, double eta) {
        return INIT_VALUE * (1 - ksi) * (1 + eta);
    }
}
