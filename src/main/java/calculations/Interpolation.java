package calculations;

import readers.DataReader;
import utils.shape_functions.ShapeFunctionDerativatives;
import utils.shape_functions.ShapeFunctions;

public class Interpolation {
    private final static int N_4 = 4; //liczba pochodnych N1/dKsi N2/dksi ... N4/dKsi

    public double[][] calculateKsiDerativatives(DataReader data) {
        int numberOfpc = (int) Math.pow(data.getNumberOfPcs(), 2);
        double[][] ksiDerativatives = new double[numberOfpc][4];

        for (int i = 0; i < numberOfpc; i++) { // idziemy po wszystkich punktach całkowania
            double id = data.getIntegrationPoints().get(i).getEta(); // 4 punkty
            ksiDerativatives[i][0] = ShapeFunctionDerativatives.ksiDerivativePC1(id); // dla 1 pkt całkowani przekazuje wartość
            ksiDerativatives[i][1] = ShapeFunctionDerativatives.ksiDerivativePC2(id);
            ksiDerativatives[i][2] = ShapeFunctionDerativatives.ksiDerivativePC3(id);
            ksiDerativatives[i][3] = ShapeFunctionDerativatives.ksiDerivativePC4(id);
        }
        return ksiDerativatives;
    }

    public double[][] calculateEtaDerativatives(DataReader data) {
        int numberOfpc = (int) Math.pow(data.getNumberOfPcs(), 2);
        double[][] etaDerativatives = new double[numberOfpc][4];

        for (int i = 0; i < numberOfpc; i++) {// idziemy po wszystkich punktach całkowania
            double id = data.getIntegrationPoints().get(i).getKsi();  // 4 punkty
            etaDerativatives[i][0] = ShapeFunctionDerativatives.etaDerivativePC1(id); // dla 1 pkt całkowani przekazuje wartość
            etaDerativatives[i][1] = ShapeFunctionDerativatives.etaDerivativePC2(id);
            etaDerativatives[i][2] = ShapeFunctionDerativatives.etaDerivativePC3(id);
            etaDerativatives[i][3] = ShapeFunctionDerativatives.etaDerivativePC4(id);
        }
        return etaDerativatives;
    }

    // liczenie funkcji kształtu
    public double[][] calculateShapeFunctionValuesMatrix(DataReader data) {
        int numberOfPcs = (int) Math.pow(data.getNumberOfPcs(), 2); //  4 f. kształtu, bo są 4 węzły
        double[][] sf = new double[numberOfPcs][numberOfPcs];
        for (int i = 0; i < numberOfPcs; i++) { // w argumentach podaje ksi i eta w zależnosci od id [-0.5773502691,0.5773502691] itd
            sf[i][0] = ShapeFunctions.PC1(data.getIntegrationPoints().get(i).getKsi(), data.getIntegrationPoints().get(i).getEta());
            // liczę funkcję kształtu podając ksi i eta do 0.25 * (1 - ksi) * (1 - eta)
            sf[i][1] = ShapeFunctions.PC2(data.getIntegrationPoints().get(i).getKsi(), data.getIntegrationPoints().get(i).getEta());
            sf[i][2] = ShapeFunctions.PC3(data.getIntegrationPoints().get(i).getKsi(), data.getIntegrationPoints().get(i).getEta());
            sf[i][3] = ShapeFunctions.PC4(data.getIntegrationPoints().get(i).getKsi(), data.getIntegrationPoints().get(i).getEta());
        }
        return sf;
    }

    public double[][] calculateShapeFunctionValuesMatrix3Points(DataReader data) {
        int numberOfPcs = (int) Math.pow(data.getNumberOfPcs3(), 2);
        double[][] sf = new double[numberOfPcs][4];

        for (int i = 0; i < numberOfPcs; i++) {
            sf[i][0] = ShapeFunctions.PC1(data.getIntegrationPoints3().get(i).getKsi(), data.getIntegrationPoints3().get(i).getEta());
            sf[i][1] = ShapeFunctions.PC2(data.getIntegrationPoints3().get(i).getKsi(), data.getIntegrationPoints3().get(i).getEta());
            sf[i][2] = ShapeFunctions.PC3(data.getIntegrationPoints3().get(i).getKsi(), data.getIntegrationPoints3().get(i).getEta());
            sf[i][3] = ShapeFunctions.PC4(data.getIntegrationPoints3().get(i).getKsi(), data.getIntegrationPoints3().get(i).getEta());
        }
        return sf;
    }

    public static double[][] calculateShapeFunctionsSurface(double[] ksi, double[] eta, DataReader data) {
        int numberOfPoints = data.getNumberOfPcs();
        int numberOfPcs = (int) Math.pow(data.getNumberOfPcs(), 2);
        double[][] N = new double[numberOfPoints][numberOfPcs]; // [2][4]

        for (int i = 0; i < 2; i++) {
            N[i][0] = ShapeFunctions.PC1(ksi[i], eta[i]); // wywolujemy dla ksi[0] i ksi[1] itd
            N[i][1] = ShapeFunctions.PC2(ksi[i], eta[i]);
            N[i][2] = ShapeFunctions.PC3(ksi[i], eta[i]);
            N[i][3] = ShapeFunctions.PC4(ksi[i], eta[i]);
        }
        return N;
    }

    //////////// dla 9 ////////////
    public double[][] calculateKsiDerativativesFor3Point(DataReader data) {
        int numberOfpcs = (int) Math.pow(data.getNumberOfPcs3(), 2);
        double[][] ksiDerativatives = new double[numberOfpcs][4]; //[9][4]

        for (int i = 0; i < numberOfpcs; i++) {
            double pcIndex = data.getIntegrationPoints3().get(i).getEta();

            ksiDerativatives[i][0] = ShapeFunctionDerativatives.ksiDerivativePC1(pcIndex);
            ksiDerativatives[i][1] = ShapeFunctionDerativatives.ksiDerivativePC2(pcIndex);
            ksiDerativatives[i][2] = ShapeFunctionDerativatives.ksiDerivativePC3(pcIndex);
            ksiDerativatives[i][3] = ShapeFunctionDerativatives.ksiDerivativePC4(pcIndex);
        }
        return ksiDerativatives;
    }

    public double[][] calculateEtaDerativativesFor3Point(DataReader data) {
        int numberOfpcs = (int) Math.pow(data.getNumberOfPcs3(), 2);
        double[][] etaDerativatives = new double[numberOfpcs][4];

        for (int i = 0; i < numberOfpcs; i++) {
            double pcIndex = data.getIntegrationPoints3().get(i).getKsi();

            etaDerativatives[i][0] = ShapeFunctionDerativatives.etaDerivativePC1(pcIndex);
            etaDerativatives[i][1] = ShapeFunctionDerativatives.etaDerivativePC2(pcIndex);
            etaDerativatives[i][2] = ShapeFunctionDerativatives.etaDerivativePC3(pcIndex);
            etaDerativatives[i][3] = ShapeFunctionDerativatives.etaDerivativePC4(pcIndex);
        }
        return etaDerativatives;
    }


}
