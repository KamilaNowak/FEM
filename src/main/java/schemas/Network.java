package schemas;

import java.util.ArrayList;
import java.util.List;

public class Network {

    private List<Node> nodes = new ArrayList<>();
    private List<Element> elements = new ArrayList<>();

    public Network(List<Node> nodes, List<Element> elements) {
        this.nodes = nodes;
        this.elements = elements;
    }

    public Network() {
    }

    public List<Node> getNodes() {
        return nodes;
    }

    public List<Element> getElements() {
        return elements;
    }


    public void printNodes() {
        for ( int i =0;i< nodes.size(); i++)
            System.out.println(i+" : (" + nodes.get(i).getX() + ", " + nodes.get(i).getY() + ")");
    }

    public void printElements() {
        for (int i =0;i< elements.size(); i++)
            System.out.print("\n[" + elements.get(i) + "]");
    }
}
