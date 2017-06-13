package sampsonLab;


import java.util.ArrayList;

/**
 * Created by dfermin on 1/6/17.
 */
public class dbNSFP_Features extends FeatureClass {

    public String mutationTaster;
    public String polyphen2_hvar;
    public String SIFT;
    public double GERP__RS;

    public dbNSFP_Features() {
        mutationTaster = "#NULL";
        polyphen2_hvar = "#NULL";
        SIFT = "#NULL";
        GERP__RS = Double.NaN;

        setDataType("MUTATIONTASTER", "VARCHAR(50)");
        setDataType("POLYPHEN2_HVAR", "VARCHAR(50)");
        setDataType("SIFT", "VARCHAR(50)");
        setDataType("GERP", "DOUBLE");
    }


    // The SIFT field is for the REF and ALT. We will record the ALT value here
    public void getSIFT(Object o) {

        if(o != null) {
            if(o.getClass().getSimpleName().equalsIgnoreCase("String")) {
                SIFT = (String) o;
            }
            else {
                ArrayList<String> inputAry = (ArrayList<String>) o;
                SIFT = inputAry.get(1);
            }
        }
    }


}
