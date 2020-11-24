package schemas;

import java.util.Arrays;
import java.util.List;

public class Element {

    private List<Integer> IDs;
    private List<Node> nodes;
    private double[][] H = new double[4][4];
    private double[][] C = new double[4][4];

    public Element() {
    }

    public List<Integer> getiDs() {
        return IDs;
    }

    public void setiDs(List<Integer> IDs) {
        this.IDs = IDs;
    }

    public List<Node> getNodes() {
        return nodes;
    }

    public void setNodes(List<Node> nodes) {
        this.nodes = nodes;
    }

    public double[][] getH() {
        return H;
    }

    public void setH(double[][] h) {
        H = h;
    }

    public double[][] getC() {
        return C;
    }

    public void setC(double[][] c) {
        C = c;
    }

    @Override
    public String toString() {
        return "Element{" +
                "iDs=" + IDs +
                ", nodes=" + nodes +
                '}';
    }
}
