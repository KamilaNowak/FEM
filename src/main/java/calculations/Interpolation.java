package calculations;

import input.DataReader;
import utils.MatrixUtils;

public class Interpolation {
    private final static int N_4 = 4; //liczba pochodnych N1/dKsi N2/dksi ... N4/dKsi
    private final static int N_9 = 9; //liczba pochodnych N1/dKsi N2/dksi ... N4/dKsi

    public double[][] calculateKsiDerativatives(DataReader data) {
        int numberOfpc = (int) Math.pow(data.getNumberOfPcs(), 2);
        double[][] ksiDerativatives = new double[numberOfpc][4];

        for (int i = 0; i < numberOfpc; i++) {
            double pcIndex = data.getIntegrationPoints().get(i).getEta();
            ksiDerativatives[i][0] = MatrixUtils.ksiDerivativePC1(pcIndex);
            ksiDerativatives[i][1] = MatrixUtils.ksiDerivativePC2(pcIndex);
            ksiDerativatives[i][2] = MatrixUtils.ksiDerivativePC3(pcIndex);
            ksiDerativatives[i][3] = MatrixUtils.ksiDerivativePC4(pcIndex);
        }
        return ksiDerativatives;
    }


    public double[][] calculateEtaDerativatives(DataReader data) {
        int numberOfpc = (int) Math.pow(data.getNumberOfPcs(), 2);
        double[][] etaDerativatives = new double[numberOfpc][4];

        for (int i = 0; i < numberOfpc; i++) {
            double pcIndex = data.getIntegrationPoints().get(i).getKsi();

            etaDerativatives[i][0] = MatrixUtils.etaDerivativePC1(pcIndex);
            etaDerativatives[i][1] = MatrixUtils.etaDerivativePC2(pcIndex);
            etaDerativatives[i][2] = MatrixUtils.etaDerivativePC3(pcIndex);
            etaDerativatives[i][3] = MatrixUtils.etaDerivativePC4(pcIndex);
        }
        return etaDerativatives;
    }


    public double[][] calculateShapeFunctionValuesMatrix(DataReader data) {
        int numberOfPcs = (int) Math.pow(data.getNumberOfPcs(), 2);
        double[][] sf = new double[numberOfPcs][4];

        for (int i = 0; i < numberOfPcs; i++) {
            sf[i][0] = MatrixUtils.shapeFunctionPC1( data.getIntegrationPoints().get(i).getKsi(), data.getIntegrationPoints().get(i).getEta());
            sf[i][1] = MatrixUtils.shapeFunctionPC2( data.getIntegrationPoints().get(i).getKsi(), data.getIntegrationPoints().get(i).getEta());
            sf[i][2] = MatrixUtils.shapeFunctionPC3( data.getIntegrationPoints().get(i).getKsi(), data.getIntegrationPoints().get(i).getEta());
            sf[i][3] = MatrixUtils.shapeFunctionPC4( data.getIntegrationPoints().get(i).getKsi(), data.getIntegrationPoints().get(i).getEta());
        }
        return sf;
    }

    //////////// dla 9 ////////////
    public double[][] calculateKsiDerativativesFor3Point(DataReader data) {
        int numberOfpcs = (int) Math.pow(data.getNumberOfPcs3(), 2);
        double[][] ksiDerativatives = new double[numberOfpcs][4]; //[9][4]

        for (int i = 0; i < numberOfpcs; i++) {
            double pcIndex = data.getIntegrationPoints3().get(i).getEta();

            ksiDerativatives[i][0] = MatrixUtils.ksiDerivativePC1(pcIndex);
            ksiDerativatives[i][1] = MatrixUtils.ksiDerivativePC2(pcIndex);
            ksiDerativatives[i][2] = MatrixUtils.ksiDerivativePC3(pcIndex);
            ksiDerativatives[i][3] = MatrixUtils.ksiDerivativePC4(pcIndex);
        }
        return ksiDerativatives;
    }

    public double[][] calculateEtaDerativativesFor3Point(DataReader data) {
        int numberOfpcs = (int) Math.pow(data.getNumberOfPcs3(), 2);
        double[][] etaDerativatives = new double[numberOfpcs][4];

        for (int i = 0; i < numberOfpcs; i++) {
            double pcIndex = data.getIntegrationPoints3().get(i).getKsi();

            etaDerativatives[i][0] = MatrixUtils.etaDerivativePC1(pcIndex);
            etaDerativatives[i][1] = MatrixUtils.etaDerivativePC2(pcIndex);
            etaDerativatives[i][2] = MatrixUtils.etaDerivativePC3(pcIndex);
            etaDerativatives[i][3] = MatrixUtils.etaDerivativePC4(pcIndex);
        }
        return etaDerativatives;
    }



}
