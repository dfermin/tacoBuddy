package sampsonLab;


import java.util.ArrayList;

/**
 * Created by dfermin on 1/6/17.
 */
public class dbNSFP_Features {

    public  String mutationTaster;
    public String polyphen2_hvar;
    public String SIFT;
    public double GERP__RS;

    public dbNSFP_Features() {
        mutationTaster = "";
        polyphen2_hvar = "";
        SIFT = "";
        GERP__RS = 0;
    }


    // The SIFT field is for the REF and ALT. We will record the ALT value here
    public void getSIFT(Object o) {
        ArrayList<String> inputAry = (ArrayList<String>) o;
        SIFT = inputAry.get(1);
    }


}
