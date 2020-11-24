package input;

import org.json.simple.JSONArray;
import org.json.simple.JSONObject;
import org.json.simple.parser.JSONParser;
import schemas.InterpolationNode;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.List;

public class DataReader {

    private double H;
    private double W;
    private double nH;
    private double nW;
    private double nE;
    private double k;
    private double c;
    private double ro;
    private double t0;

    private List<InterpolationNode> integrationPoints = new ArrayList<>();
    private List<InterpolationNode> integrationPoints3 = new ArrayList<>();
    private List<Double> pc2 = new ArrayList<>();
    private int numberOfPcs = 0;

    private List<Double> pc3 = new ArrayList<>();
    private int numberOfPcs3 = 0;

    private List<Double> weights = new ArrayList<>();
    private int numberOfWeights = 0;

    private List<Double> weights3 = new ArrayList<>();
    private int numberOfWeights3 = 0;

    public DataReader() throws FileNotFoundException {
        this.readConstantData();
    }

    private void readConstantData() throws FileNotFoundException {
        JSONParser parser = new JSONParser();

        Object input;
        try {
            input = parser.parse(new FileReader(System.getProperty("user.dir") + "/src/main/java/input/data.json"));
        } catch (Exception e) {
            throw new FileNotFoundException();
        }
        JSONObject constants = (JSONObject) input;
        this.H = (double) constants.get("H");
        this.W = (double) constants.get("W");
        this.nH = (double) constants.get("nH");
        this.nW = (double) constants.get("nW");
        this.k = (double) constants.get("k");
        this.nE = (double) constants.get("nE");
        this.c = (double) constants.get("c");
        this.ro = (double) constants.get("ro");
        this.t0 = (double) constants.get("t0");

        JSONArray pcInput = (JSONArray) constants.get("pc2");
        this.numberOfPcs = pcInput.size();
        for (Object value : pcInput) this.pc2.add((Double) value);

        JSONArray pcInput3 = (JSONArray) constants.get("pc3");
        this.numberOfPcs3 = pcInput3.size();
        for (Object value : pcInput3) this.pc3.add(Double.parseDouble(value.toString()));

        JSONArray weightsInput = (JSONArray) constants.get("weights");
        this.numberOfWeights = weightsInput.size();
        for (Object value : weightsInput) this.weights.add((Double) value);

        JSONArray weightsInput3 = (JSONArray) constants.get("weights3");
        this.numberOfWeights3 = weightsInput3.size();
        for (Object value : weightsInput3) this.weights3.add((Double) value);

        for (int i = 0; i < numberOfPcs; i++) {
            for (int j = 0; j < numberOfPcs; j++)
                integrationPoints.add(new InterpolationNode(this.pc2.get(j), this.pc2.get(i)));
        }

        for (int i = 0; i < numberOfPcs3; i++) {
            for (int j = 0; j < numberOfPcs3; j++) {
                InterpolationNode in = new InterpolationNode(this.pc3.get(j), this.pc3.get(i));
                if (this.pc3.get(j) == 0 || this.pc3.get(i) == 0) {
                    in.setWeight(getWeights3().get(1));
                } else {
                    in.setWeight(getWeights3().get(0));
                }
                integrationPoints3.add(in);
            }
        }
        System.out.println("getWeights3() : " + getWeights3());

    }

    public List<Double> getWeights3() {
        return weights3;
    }

    public double getnE() {
        return nE;
    }

    public double getH() {
        return H;
    }

    public double getW() {
        return W;
    }

    public double getnH() {
        return nH;
    }

    public double getnW() {
        return nW;
    }

    public double getK() {
        return k;
    }

    public List<Double> getPc3() {
        return pc3;
    }

    public int getNumberOfPcs3() {
        return numberOfPcs3;
    }

    public List<InterpolationNode> getIntegrationPoints() {
        return integrationPoints;
    }


    public List<InterpolationNode> getIntegrationPoints3() {
        return integrationPoints3;
    }

    public void setIntegrationPoints3(List<InterpolationNode> integrationPoints3) {
        this.integrationPoints3 = integrationPoints3;
    }

    public List<Double> getPc2() {
        return pc2;
    }

    public void setPc2(List<Double> pc) {
        this.pc2 = pc;
    }

    public int getNumberOfPcs() {
        return numberOfPcs;
    }

    public void setNumberOfPcs(int numberOfPcs) {
        this.numberOfPcs = numberOfPcs;
    }

    public List<Double> getWeights() {
        return weights;
    }

    public void setWeights(List<Double> weights) {
        this.weights = weights;
    }

    public int getNumberOfWeights() {
        return numberOfWeights;
    }

    public void setNumberOfWeights(int numberOfWeights) {
        this.numberOfWeights = numberOfWeights;
    }

    public double getC() {
        return c;
    }

    public double getRo() {
        return ro;
    }

    public double getT0() {
        return t0;
    }

    @Override
    public String toString() {
        return "DataReader{" +
                "\nH=" + H +
                "\n W=" + W +
                "\nnH=" + nH +
                "\nnW=" + nW +
                "\nnE=" + nE +
                "\n k=" + k +
                "\n c=" + c +
                "\n ro=" + ro +
                "\n t0=" + t0 +
                "\nintegrationPoints=" + integrationPoints +
                "\nintegrationPoints3=" + integrationPoints3 +
                "\npc2=" + pc2 +
                "\nnumberOfPcs=" + numberOfPcs +
                "\npc3=" + pc3 +
                "\nnumberOfPcs3=" + numberOfPcs3 +
                "\nweights=" + weights +
                "\n numberOfWeights=" + numberOfWeights +
                "\neights3=" + weights3 +
                "\nnumberOfWeights3=" + numberOfWeights3 +
                '}';
    }
}
