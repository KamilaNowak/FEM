package schemas;

public class Node {
    private int id;
    private double x;
    private double y;
    private double t;

    public Node(double x, double y) {
        this.x = x;
        this.y = y;
    }

    public Node() {
    }

    public void setId(int id) {
        this.id = id;
    }

    public double getX() {
        return x;
    }

    public double getY() {
        return y;
    }

    public double getT() {
        return t;
    }

    public void setT(double t) {
        this.t = t;
    }

    @Override
    public String toString() {
        return "Node{" +
                "id=" + id +
                ", x=" + x +
                ", y=" + y +
                ", t=" + t +
                '}';
    }
}
