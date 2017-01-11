package sampsonLab;

import com.google.common.base.Joiner;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.jexl3.*;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

/**
 * Created by dfermin on 1/5/17.
 */
public class VariantInfo {
    public String chr;
    public int pos;
    public String REF;
    public String ALT;
    public String snp_id;
    public int passedFilter; // 0 = false, 1 = true

    public double svmProb;

    // These object only get initialized if the user decides to use them
    public ESP_Features ESP = null;
    public EFF_Features EFF = null;
    public dbNSFP_Features dbNSFP = null;



    public VariantInfo(String chr, int pos, String ref, String alt) {
        this.chr = chr;
        this.pos = pos;
        this.REF = ref;
        this.ALT = alt;
        this.svmProb = 0.0;
        this.passedFilter = 0;

        // construct snp_id string
        snp_id = chr + ":" + String.valueOf(pos) + ":" + REF + ">" + ALT;

        ESP = new ESP_Features();
        dbNSFP = new dbNSFP_Features();
    }

    public String getID() { return snp_id; }
    public void setSvmProb(double d) { this.svmProb = d; }




    // Function records the requested data from the attributes of the variant object
    public boolean fetchFeature(String feat, VariantContext vc) {

        if(feat.startsWith("ESP_")) {
            ESP = new ESP_Features();

            String aa = vc.getAttributeAsString("ESP_AA_AC", ".").replaceAll("\\[", "").replaceAll("\\]", "");
            String ea = vc.getAttributeAsString("ESP_EA_AC", ".").replaceAll("\\[", "").replaceAll("\\]", "");
            ESP.ESP_AA_AC = ESP.calcMAF(aa);
            ESP.ESP_EA_AC = ESP.calcMAF(ea);
            ESP.ESP_MAX_AA_EA = Math.max(ESP.ESP_AA_AC, ESP.ESP_EA_AC);

            ESP.setESP_MAF(vc.getAttributeAsString("ESP_MAF", ".").replaceAll("\\[", "").replaceAll("\\]", ""));
        }

        // Search for this feature in the VariantInfoNameMap
        // If you get a match, record the relevant data for that feature
        if( tacoBuddy.VCF_Info_Name_Map.inMap(feat.toUpperCase()) ) {
            String k = tacoBuddy.VCF_Info_Name_Map.getValue(feat.toUpperCase());
            //Object v = vc.getAttribute(k);

            String v_str = vc.getAttributeAsString(k, "#NULL");
            v_str = v_str.replaceAll("\\[", "").replaceAll("\\]", "");

            if( !v_str.equalsIgnoreCase("#NULL") ) {

                if (feat.equalsIgnoreCase("EFF")) {
                    EFF = new EFF_Features(v_str);
                }
                if (feat.equalsIgnoreCase("MutationTaster")) {
                    if (dbNSFP == null) dbNSFP = new dbNSFP_Features();
                    dbNSFP.mutationTaster = v_str;
                }

                if (feat.equalsIgnoreCase("GERP")) {
                    if (dbNSFP == null) dbNSFP = new dbNSFP_Features();
                    dbNSFP.GERP__RS = Double.valueOf(v_str);
                }

                if (feat.equalsIgnoreCase("Polyphen2_Hvar")) {
                    if (dbNSFP == null) dbNSFP = new dbNSFP_Features();
                    dbNSFP.polyphen2_hvar = v_str;
                }

                if (feat.equalsIgnoreCase("sift")) {
                    if (dbNSFP == null) dbNSFP = new dbNSFP_Features();
                    dbNSFP.SIFT = v_str;
                }
            }
        }

        return true; // if you got here then you had no problems in this function
    }


    // Function applies the user's filter to this variant. Calling this function assumes VariantInfo::fetchFeature()
    // ran successfully. The curTS_ID is only needed if the user requested to evaluate the 'EFF' info field.
    public void filter(String jexl_filter_str, String curTS_ID) {

        // The ++ string breaks jexl so we strip it out. It should only occur in the GERP++ name anyways
        if(jexl_filter_str.contains("++")) {
            String tmp = jexl_filter_str.replaceAll("\\+\\+", "");
            jexl_filter_str = tmp;
            tmp = null;
        }

        // The '/' string breaks jexl so we strip it out. It should only occur in the GERP++ name anyways
        if(jexl_filter_str.contains("/")) {
            String tmp = jexl_filter_str.replaceAll("/", "");
            jexl_filter_str = tmp;
            tmp = null;
        }

        jexl_filter_str = jexl_filter_str.replaceAll(" and ", " && ");
        jexl_filter_str = jexl_filter_str.replaceAll(" or ", " || ");


        JexlEngine jexl = new JexlBuilder().create(); // Create a jexl engine
        JexlExpression expr = jexl.createExpression(jexl_filter_str); // define the expression you want to test/use

        // Create a MapContext object and populate it with the variables that are in VariantInfo objects
        JexlContext jc = new MapContext();

        if(dbNSFP != null) prepContextMap(dbNSFP, jc, curTS_ID);
        if(ESP != null) prepContextMap(ESP, jc, curTS_ID);
        if(EFF != null) prepContextMap(EFF, jc, curTS_ID);

        boolean retVal = false;
        retVal = (Boolean) expr.evaluate(jc);

        if(retVal) passedFilter = 1;
    }



    // Function to populate the jexlContext map with the variables the user has elected to filter on
    private void prepContextMap(Object o, JexlContext jc, String curTS_ID) {

        String dataType = o.getClass().getSimpleName();

        if(dataType.equalsIgnoreCase("dbNSFP_Features")) {
            jc.set("GERP", dbNSFP.GERP__RS);
            jc.set("MUTATIONTASTER", dbNSFP.mutationTaster);
            jc.set("POLYPHEN2_HVAR", dbNSFP.polyphen2_hvar);
            jc.set("SIFT", dbNSFP.SIFT);
        }

        if(dataType.equalsIgnoreCase("ESP_Features")) {
            jc.set("ESP_MAX_AA_EA", ESP.ESP_MAX_AA_EA);
            jc.set("ESP_MAF", ESP.ESP_MAF);
        }

        if(dataType.equalsIgnoreCase("EFF_Features")) {
            String jj = EFF.findTS(curTS_ID);
            jc.set("EFF", EFF.findTS(curTS_ID));
        }
    }

    // Function returns the requested features for this variant.
    public String returnSummaryString(String geneId, boolean printThisVariant) {
        String ret = "";

        if(printThisVariant == false) return "";

        HashMap<String, String> transcriptMap = null;
        HashMap<String, String> returnValueMap = new HashMap<String, String>();

        for(String feat : globalFunctions.featureSet) {

            returnValueMap.put(feat, "#NULL"); // Initialize map and update it as data is encountered

            // Special case where you have to iterate over the transcripts.
            if (feat.equalsIgnoreCase("EFF")) {
                transcriptMap = EFF.eff;
            }
            if (feat.equalsIgnoreCase("MUTATIONTASTER")) {
                returnValueMap.put(feat, dbNSFP.mutationTaster);
                continue;
            }
            if (feat.equalsIgnoreCase("POLYPHEN2_HVAR")) {
                returnValueMap.put(feat, dbNSFP.polyphen2_hvar);
                continue;
            }
            if (feat.equalsIgnoreCase("SIFT")) {
                returnValueMap.put(feat, dbNSFP.SIFT);
                continue;
            }
            if (feat.equalsIgnoreCase("GERP")) {
                returnValueMap.put(feat, dbl2str(dbNSFP.GERP__RS, 4));
                continue;
            }

            if(feat.equalsIgnoreCase("ESP_AA_AC")) {
                returnValueMap.put(feat, dbl2str(ESP.ESP_AA_AC, 4) );
                continue;
            }
            if(feat.equalsIgnoreCase("ESP_EA_AC")) {
                returnValueMap.put(feat, dbl2str(ESP.ESP_EA_AC, 4) );
                continue;
            }
            if(feat.equalsIgnoreCase("ESP_MAF")) {
                returnValueMap.put(feat, dbl2str(ESP.ESP_MAF, 4) );
                continue;
            }
            if(feat.equalsIgnoreCase("ESP_MAX_AA_EA")) {
                returnValueMap.put(feat, dbl2str(ESP.ESP_MAX_AA_EA, 4) );
                continue;
            }
        }

        // Now construct the return string based upon the data in returnValueMap.
        ArrayList<String> dataList = new ArrayList<String>();
        for(String feat : globalFunctions.featureSet) {
            if(feat.equalsIgnoreCase("EFF")) continue;
            String tmp = feat + "=" + returnValueMap.get(feat);
            dataList.add(tmp);
        }
        String data_line = Joiner.on("\t").join(dataList) + "\n";


        if(transcriptMap != null) {
            // Special case, need to report each transcript row by row
            String base = Integer.toString(this.passedFilter) + "\t" + snp_id + "\t" + geneId + "\t";
            for(String curTS : transcriptMap.keySet()) {
                ret += base + curTS + "\tEFF=" + transcriptMap.get(curTS) + "\t" + data_line;
            }
        }
        else {
            ret = Integer.toString(this.passedFilter) + "\t" + snp_id + "\t" + geneId + "\t" + data_line;
        }

        return ret;
    }



    // Manually round the given double to a specific number of digits and return a String with it's value
    private static String dbl2str(double d, int precision) {
        String t = "%." + Integer.toString(precision) + "f";
        return String.format(t, d);
    }
}
