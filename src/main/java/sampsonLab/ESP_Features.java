package sampsonLab;

import java.util.ArrayList;
import java.util.Collections;


/**
 * Created by dfermin on 1/5/17.
 */
public class ESP_Features extends FeatureClass {

    public double ESP_AA_AC, ESP_EA_AC, ESP_MAF;
    public double ESP_MAX_AA_EA;

    public ESP_Features() {
        ESP_AA_AC = -1;
        ESP_EA_AC = -1;
        ESP_MAF = -1;
        ESP_MAX_AA_EA = -1;
    }


    public double calcMAF(String AC_str) {
        ArrayList<Double> counts = new ArrayList<Double>();
        ArrayList<Double> AF = new ArrayList<Double>();
        double sum = 0;

        if(AC_str.equalsIgnoreCase(".")) return 0.0;

        // record all of the allele counts
        for(String s : AC_str.split(",")) {
            double d = Double.valueOf(s);
            counts.add(d);
            sum += d;
        }

        // Doing stuff in log-scale to avoid buffer underflow error
        for(double d : counts) AF.add( Math.log(d / sum) );

        // Sort AF from high to low and accept the second element as the minor allele frequency (MAF)
        Collections.sort(AF);
        Collections.reverse(AF);

        double MAF_ret = 0;
        if(AF.size() > 1) MAF_ret = Math.exp(AF.get(1)); // convert back to linear-scale!!!
        else MAF_ret = Math.exp(AF.get(0));  // convert back to linear-scale!!!

        return MAF_ret;
    }

    // Parses the ESP_MAF field from the given string
    public void setESP_MAF(String AC_str) {
        if( !AC_str.equalsIgnoreCase(".") ) {
            ArrayList<Double> AF = new ArrayList<Double>();
            for(String s : AC_str.split(",")) {
                double d = Double.valueOf(s);
                AF.add(d);
            }
            Collections.sort(AF);
            Collections.reverse(AF);
            double tmp = 0;
            if(AF.size() > 1) this.ESP_MAF = AF.get(1);
            else this.ESP_MAF = AF.get(0);
        }
    }

}
