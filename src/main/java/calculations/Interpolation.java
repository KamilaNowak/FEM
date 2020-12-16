package calculations;

import input.DataReader;
import utils.Derativatives;
import utils.ShapeFunctions;

public class Interpolation {
    private final static int N_4 = 4; //liczba pochodnych N1/dKsi N2/dksi ... N4/dKsi
    private final static int N_9 = 9; //liczba pochodnych N1/dKsi N2/dksi ... N4/dKsi

    public double[][] calculateKsiDerativatives(DataReader data) {
        int numberOfpc = (int) Math.pow(data.getNumberOfPcs(), 2);
        double[][] ksiDerativatives = new double[numberOfpc][4];

        for (int i = 0; i < numberOfpc; i++) {
            double id = data.getIntegrationPoints().get(i).getEta();
            ksiDerativatives[i][0] = Derativatives.ksiDerivativePC1(id);
            ksiDerativatives[i][1] = Derativatives.ksiDerivativePC2(id);
            ksiDerativatives[i][2] = Derativatives.ksiDerivativePC3(id);
            ksiDerativatives[i][3] = Derativatives.ksiDerivativePC4(id);
        }
        return ksiDerativatives;
    }


    public double[][] calculateEtaDerativatives(DataReader data) {
        int numberOfpc = (int) Math.pow(data.getNumberOfPcs(), 2);
        double[][] etaDerativatives = new double[numberOfpc][4];

        for (int i = 0; i < numberOfpc; i++) {
            double id = data.getIntegrationPoints().get(i).getKsi();

            etaDerativatives[i][0] = Derativatives.etaDerivativePC1(id);
            etaDerativatives[i][1] = Derativatives.etaDerivativePC2(id);
            etaDerativatives[i][2] = Derativatives.etaDerivativePC3(id);
            etaDerativatives[i][3] = Derativatives.etaDerivativePC4(id);
        }
        return etaDerativatives;
    }


    public double[][] calculateShapeFunctionValuesMatrix(DataReader data) {
        int numberOfPcs = (int) Math.pow(data.getNumberOfPcs(), 2);
        double[][] sf = new double[numberOfPcs][4];
        for (int i = 0; i < numberOfPcs; i++) {
            sf[i][0] = ShapeFunctions.PC1(data.getIntegrationPoints().get(i).getKsi(), data.getIntegrationPoints().get(i).getEta());
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

    public static double[][] calculateShapeFunctionsSurface(double[] ksi, double[] eta) {
        double[][] N = new double[2][4];

        for (int i = 0; i < 2; i++) {
            N[i][0] = ShapeFunctions.PC1(ksi[i], eta[i]);
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

            ksiDerativatives[i][0] = Derativatives.ksiDerivativePC1(pcIndex);
            ksiDerativatives[i][1] = Derativatives.ksiDerivativePC2(pcIndex);
            ksiDerativatives[i][2] = Derativatives.ksiDerivativePC3(pcIndex);
            ksiDerativatives[i][3] = Derativatives.ksiDerivativePC4(pcIndex);
        }
        return ksiDerativatives;
    }

    public double[][] calculateEtaDerativativesFor3Point(DataReader data) {
        int numberOfpcs = (int) Math.pow(data.getNumberOfPcs3(), 2);
        double[][] etaDerativatives = new double[numberOfpcs][4];

        for (int i = 0; i < numberOfpcs; i++) {
            double pcIndex = data.getIntegrationPoints3().get(i).getKsi();

            etaDerativatives[i][0] = Derativatives.etaDerivativePC1(pcIndex);
            etaDerativatives[i][1] = Derativatives.etaDerivativePC2(pcIndex);
            etaDerativatives[i][2] = Derativatives.etaDerivativePC3(pcIndex);
            etaDerativatives[i][3] = Derativatives.etaDerivativePC4(pcIndex);
        }
        return etaDerativatives;
    }


}
