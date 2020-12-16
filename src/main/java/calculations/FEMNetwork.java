package calculations;

import input.DataReader;
import schemas.Element;
import schemas.Network;
import schemas.Node;

import java.util.ArrayList;
import java.util.List;

public class FEMNetwork {

    private List<Node> nodes;
    private List<Element> elements;

    public void createNodes(DataReader dataReader, List<Node> resultNodes) {

        double H = dataReader.getH();
        double W = dataReader.getW();
        double nW = dataReader.getNW();
        double nH = dataReader.getNH();

        double dX = W / (nW - 1);
        double dY = H / (nH - 1);
        int id = 0;

        Node node;
        for (double i = 0; i < nW; i++) {
            for (double j = 0; j < nH; j++) {

                id++;
                double x = i * dX;
                double y = j * dY;

                node = new Node(x, y);
                node.setId(id);
                node.setT(dataReader.getT0());

                //dla Hbc
                if (x == 0 || x == W || y == 0 || y == H )
                    node.setBc(true);

                resultNodes.add(node);
            }
        }
    }

    public void createElements(DataReader dataReader, List<Element> resultElements) {
        double nH = dataReader.getNH();
        double nE = dataReader.getNE();
        int elem = 0;

        for (int i = 1; i <= nE + 2; i++) {
            if (i % nH == 0) elem++;
        }

        for (int i = 1; i <= nE + 2; i++) {
            if (i % nH == 0) continue;

            Element element = new Element();
            List<Integer> elementIDs = new ArrayList<>();
            List<Node> nodes = new ArrayList<>();

            int id1 = i;
            int id2 = (int) (i + nH);
            int id3 = (int) (i + nH + 1);
            int id4 = i + 1;

            elementIDs.add(id1);
            elementIDs.add(id2);
            elementIDs.add(id3);
            elementIDs.add(id4);

            nodes.add(this.nodes.get(id1 - 1));
            nodes.add(this.nodes.get(id2 - 1));
            nodes.add(this.nodes.get(id3 - 1));
            nodes.add(this.nodes.get(id4 - 1));

            element.setIDs(elementIDs);
            element.setNodes(nodes);
            resultElements.add(element);
        }
    }

    public Network createNetwork(DataReader data) {
        Network network = new Network();

        nodes = network.getNodes();
        elements = network.getElements();

        this.createNodes(data, nodes);
        this.createElements(data, elements);
        return network;
    }
}
