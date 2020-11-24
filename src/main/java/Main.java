import calculations.FEMNetwork;
import calculations.Interpolation;
import calculations.MatrixCalculates;
import input.DataReader;
import printable.Print;
import schemas.Element;
import schemas.Network;
import sun.nio.ch.Net;
import utils.MatrixUtils;

import java.io.FileNotFoundException;

public class Main {
    public static void main(String[] args) throws FileNotFoundException {

        MatrixCalculates matrixCalculates = new MatrixCalculates();

        DataReader data = new DataReader();
        System.out.println("Input Data: ");
        System.out.println(data);

        FEMNetwork femNetwork = new FEMNetwork();

        Network network = femNetwork.createNetwork(data);
        network.printNodes();
        network.printElements();

        Interpolation interpolation = new Interpolation();

        double[][] ksi = interpolation.calculateKsiDerativatives(data);
        double[][] eta = interpolation.calculateEtaDerativatives(data);

        System.out.println("\n\nPochodne po ksi:");
        Print.printMatrixNxM(ksi, ksi[0].length, ksi[1].length);
        System.out.println("\nPochodne po eta:");
        Print.printMatrixNxM(eta, eta[0].length, eta[1].length);

        System.out.println("\n\nFunkcja kształtu:");
        double[][] shapeFunction = interpolation.calculateShapeFunctionValuesMatrix(data);
        Print.printMatrixNxM(shapeFunction, shapeFunction[0].length, shapeFunction[1].length);

        System.out.println("\nMacierze H dla każdego elementu 2p:");
        for (Element element : network.getElements()) {
            double[][] H = matrixCalculates.calculateHMatrix(element, ksi, eta, data);
            System.out.println("\n H dla " + element.getNodes());
            Print.printMatrixNxM(H, 4, 4);
        }
        System.out.println("\nMacierze C dla każdego elementu 2p:");
        for (Element element : network.getElements()) {
            double[][] C = matrixCalculates.generateC(element, shapeFunction, ksi, eta, data);
            System.out.println("\n C dla " + element.getNodes());
            Print.printMatrixNxM(C, 4, 4);
        }

        System.out.println("\n Global H Matrix: ");
        double[][] H_global = matrixCalculates.calculateGlobalHMatrix(network);
        Print.printMatrixNxM(H_global, H_global.length, H_global.length);

        System.out.println("\n Global C Matrix: ");
        double[][] C_global = matrixCalculates.calculateGlobalCMatrix(network);
        Print.printMatrixNxM(C_global, C_global.length, C_global.length);

        System.out.println("\n\n----------------------------------------------------- 3 p ---------------------------------------------------------------\n\n");
        //////// dla 3p ////////
        System.out.println("KSi pochodne dla 3p");
        double[][] rKsi = interpolation.calculateKsiDerativativesFor3Point(data);
        Print.printMatrixNxM(rKsi, 9, 4);

        System.out.println("Eta pochodne dla 3p");
        double[][] rEta = interpolation.calculateEtaDerativativesFor3Point(data);
        Print.printMatrixNxM(rEta, 9, 4);

        System.out.println("\nMacierze H dla każdego elementu 3p:");
        for (Element element : network.getElements()) {
            double[][] H3 = matrixCalculates.calculateHMatrix3Points(element, rKsi, rEta, data);
            System.out.println("\n H dla " + element.getNodes());
            Print.printMatrixNxM(H3, 4, 4);
        }
        System.out.println("\n Global H Matrix: ");
        double[][] H3_global = matrixCalculates.calculateGlobalHMatrix(network);
        Print.printMatrixNxM(H3_global, H3_global.length, H3_global.length);

        System.out.println("\nMacierze C dla każdego elementu 3p:");
        for (Element element : network.getElements()) {
            double[][] C3 = matrixCalculates.generateCFor3Points(element, shapeFunction, rKsi, rEta, data);
            System.out.println("\n C dla " + element.getNodes());
            Print.printMatrixNxM(C3, 4, 4);
        }
        System.out.println("\n Global C Matrix 3p: ");
        double[][] C_global3 = matrixCalculates.calculateGlobalCMatrix(network);
        Print.printMatrixNxM(C_global3, C_global3.length, C_global3.length);


    }
}
