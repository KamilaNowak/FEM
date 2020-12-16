package input;

import lombok.*;
import org.json.simple.JSONArray;
import org.json.simple.JSONObject;
import org.json.simple.parser.JSONParser;
import schemas.InterpolationNode;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.List;

@Getter
@Setter
@AllArgsConstructor
@ToString
@Data
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
    private double alfa;
    private double tempAlfa;
    private double alfaForP;
    private double simulationTime;
    private double simulationStepTime;

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

    public DataReader() throws FileNotFoundException { this.readConstantData(); }

    private void readConstantData() throws FileNotFoundException {
        JSONParser parser = new JSONParser();
        Object input;

        try { input = parser.parse(new FileReader(System.getProperty("user.dir") + "/src/main/java/input/data.json")); }
        catch (Exception e) { throw new FileNotFoundException(); }

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
        this.alfa = (double) constants.get("alfa");
        this.tempAlfa = (double) constants.get("t_alfa");
        this.alfaForP = (double) constants.get("alfa_for_P");
        this.simulationTime = (double) constants.get("simTime");
        this.simulationStepTime= (double) constants.get("simStepTime");

        JSONArray pcInput = (JSONArray) constants.get("pc2");
        JSONArray pcInput3 = (JSONArray) constants.get("pc3");
        JSONArray weightsInput = (JSONArray) constants.get("weights");
        JSONArray weightsInput3 = (JSONArray) constants.get("weights3");

        this.numberOfPcs = pcInput.size();
        this.numberOfPcs3 = pcInput3.size();
        this.numberOfWeights = weightsInput.size();
        this.numberOfWeights3 = weightsInput3.size();

        for (Object value : pcInput) this.pc2.add((Double) value);
        for (Object value : pcInput3) this.pc3.add(Double.parseDouble(value.toString()));
        for (Object value : weightsInput) this.weights.add((Double) value);
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
      //  System.out.println("getWeights3() : " + getWeights3());

    }
}
