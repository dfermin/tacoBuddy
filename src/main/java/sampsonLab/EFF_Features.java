package sampsonLab;

import com.google.common.base.Joiner;
import com.google.common.collect.ImmutableMap;
import jdk.nashorn.internal.ir.annotations.Immutable;

import java.util.*;

/**
 * Created by dfermin on 1/6/17.
 */
public class EFF_Features extends FeatureClass {

    public HashMap<String, String> eff; // k = transcript ID, v = status_for_this_transcript
    public HashMap<String, Integer> impactMap; // k = transcript ID, v = the impact of this variant (HIGH=3, MODERATE=2, LOW=1)
    public HashMap<String, String> protChange; // k = transcript ID, v = amino acid change (if any)

    static final Map<String, String> aaMap = ImmutableMap.<String, String>builder()
            .put("Ala", "A")
            .put("Arg", "R")
            .put("Asn", "N")
            .put("Asp", "D")
            .put("Cys", "C")
            .put("Gln", "Q")
            .put("Glu", "E")
            .put("Gly", "G")
            .put("His", "H")
            .put("Ile", "I")
            .put("Leu", "L")
            .put("Lys", "K")
            .put("Met", "M")
            .put("Phe", "F")
            .put("Pro", "P")
            .put("Ser", "S")
            .put("Thr", "T")
            .put("Trp", "W")
            .put("Tyr", "Y")
            .put("Val", "V")
            .build();

    public EFF_Features(String ss) {

        eff = new HashMap<String, String>();
        protChange = new HashMap<String, String>();
        impactMap = new HashMap<String, Integer>();

        // this array will contain 1 entry per transcript affected by the variant call
        String ss2 = ss.replaceAll("[\\[\\]]+", "");
        String[] inputAry = ss2.split(",");

        for(String s : inputAry) {

            String[] tmpAry = s.split("\\|");

            String label = tmpAry[0].substring(0,tmpAry[0].indexOf('(')).trim(); // gets the "assigned effect" of this variant on the given gene.
            String tmp2 = label.replaceAll("_variant", "");
            label = tmp2;

            String transId = "";
            String protAA = "#NULL";
            String impactValue = "#NULL";
            for(int j = 0; j < tmpAry.length; j++) {

                if( j == 0 ) { //should be the impact of this variant according to http://snpeff.sourceforge.net/VCFannotationformat_v1.0.pdf
                    impactValue = tmpAry[j].split("\\(")[1];
                    continue;
                }

                if( tmpAry[j].startsWith("ENST00") ) {
                    transId = tmpAry[j].replaceAll("\\.\\d+$", "");
                }
                if( tmpAry[j].startsWith("p.")) { // amino acid level change information
                    protAA = parseProteinChange(tmpAry[j]);
                }
            }

            if(transId.startsWith("ENST00")) {
                eff.put(transId, label);
                protChange.put(transId, protAA);

                if(impactValue.equalsIgnoreCase("HIGH")) impactMap.put(transId, 3);
                else if(impactValue.equalsIgnoreCase("MODERATE")) impactMap.put(transId, 2);
                else if(impactValue.equalsIgnoreCase("LOW")) impactMap.put(transId, 1);
                else impactMap.put(transId, 0);
            }
        }
    }


    // Function returns a string indicating the change the variant causes at the protein level
    private String parseProteinChange(String ss) {
        String ret = ss.split("/")[0].substring(2);

        for(String aa3 : this.aaMap.keySet()) {
            if(ret.contains(aa3)) {
                String x1 = ret.replaceAll(aa3, this.aaMap.get(aa3));
                ret = x1;
            }
        }

        if(ret.endsWith("fs")) { // frame shift
            ret = "-" + ret.replaceAll("fs", "");
        }

        return ret;
    }


    // Function returns TRUE if the impact value for the given transcript matches the value passed to this function
    public boolean checkEFFimpact(String curTS, String cutOff) {
        boolean ret = false;

        int cutOffNum = 0;
        if(cutOff.equalsIgnoreCase("HIGH")) cutOffNum = 3;
        if(cutOff.equalsIgnoreCase("MODERATE")) cutOffNum = 2;
        if(cutOff.equalsIgnoreCase("LOW")) cutOffNum = 1;

        if(impactMap.containsKey(curTS)){
            if(impactMap.get(curTS) >= cutOffNum) ret = true;
        }

        return ret;
    }


    // Function returns TRUE if the data in EFF_Features says this variant has no affect on the phenotype
    public boolean isSynonymousVariant(String ts) {
        boolean ret = true;

        if(eff.containsKey(ts)) {
            if( !eff.get(ts).equalsIgnoreCase("synonymous") ) ret = false;
        }
        return(ret);
    }


    // Function converts the impact numbers to string terms
    private String impactNum2str(String ts) {
        String ret = ".";

        if(impactMap.containsKey(ts)) {
            int I = impactMap.get(ts);
            if(I == 1) ret = "LOW";
            if(I == 2) ret = "MODERATE";
            if(I == 3) ret = "HIGH";
        }
        return ret;
    }


    // Function returns the value of 'eff' map for the given key 'k'.
    // If no match is found, null is returned
    public String findTS(String search_str) {
        String ret = ".";

        if(eff.containsKey(search_str.replaceAll("\\.\\d+$", ""))) {
            ret = eff.get(search_str);
            ret += ":" + protChange.get(search_str) + ":" + impactNum2str(search_str);
        }
        return ret;
    }


    // Function iterates over the values in 'eff' and concatenates the unique ones into an array
    // that jexl can interpret. This function is for gene-level output.
    public String returnJexlArray() {
        String ret = "";
        SortedSet<String> SS = new TreeSet<String>();

        for(String v : eff.values()) {

            for(String i : v.split(",")) {
                i = "'" + i.trim() + "'";
                SS.add(i);
            }
        }

        ret = "[" + Joiner.on(",").join(SS) + "]";
        return ret;
    }
}
