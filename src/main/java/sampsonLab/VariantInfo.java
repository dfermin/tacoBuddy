package sampsonLab;

import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.jexl3.*;

import java.util.ArrayList;

/**
 * Created by dfermin on 1/5/17.
 */
public class VariantInfo {
    public static String chr;
    public static int pos;
    public static String REF;
    public static String ALT;
    public static String snp_id;


    // These object only get initialized if the user decides to use them
    static public ESP_Features ESP = null;
    static public EFF_Features EFF = null;
    static public dbNSFP_Features dbNSFP = null;

    public VariantInfo(String chr, int pos, String ref, String alt) {
        this.chr = chr;
        this.pos = pos;
        this.REF = ref;
        this.ALT = alt;

        // construct snp_id string
        snp_id = chr + ":" + String.valueOf(pos) + ":" + REF + ">" + ALT;
    }

    public String getID() { return snp_id; }



    // Function records the requested data from the attributes of the variant object
    public boolean fetchFeature(String feat, VariantContext vc) {

        if(feat.startsWith("ESP_")) {
            ESP = new ESP_Features();

            String aa = vc.getAttributeAsString("ESP_AA_AC", ".").replaceAll("\\[", "").replaceAll("\\]", "");
            String ea = vc.getAttributeAsString("ESP_EA_AC", ".").replaceAll("\\[", "").replaceAll("\\]", "");
            ESP.ESP_AA_AC = ESP.calcMAF(aa);
            ESP.ESP_EA_AC = ESP.calcMAF(ea);

            // This is a custom calculation you have to perform
            if(globalFunctions.customCalcsMap.containsKey(feat)) {

                if (feat.equalsIgnoreCase("ESP_MAX_AA_EA")) {
                    ESP.ESP_MAX_AA_AE = Math.max(ESP.ESP_AA_AC, ESP.ESP_EA_AC);
                }
            }
        }

        // Search for this feature in the VariantInfoNameMap
        // If you get a match, record the relevant data for that feature
        if( tacoBuddy.VCF_Info_Name_Map.inMap(feat.toUpperCase()) ) {
            String k = tacoBuddy.VCF_Info_Name_Map.getValue(feat.toUpperCase());
            Object v = vc.getAttribute(k);

            if(v != null) {
                if (feat.equalsIgnoreCase("EFF")) {
                    EFF = new EFF_Features(v);
                }
                if (feat.equalsIgnoreCase("MutationTaster")) {
                    if (dbNSFP == null) dbNSFP = new dbNSFP_Features();
                    dbNSFP.mutationTaster = (String) v;
                }

                if (feat.equalsIgnoreCase("GERP++")) {
                    if (dbNSFP == null) dbNSFP = new dbNSFP_Features();
                    dbNSFP.GERP__RS = Double.valueOf((String) v);
                }

                if (feat.equalsIgnoreCase("Polyphen2_Hvar")) {
                    if (dbNSFP == null) dbNSFP = new dbNSFP_Features();
                    dbNSFP.polyphen2_hvar = ((ArrayList<String>) v).get(0);
                }

                if (feat.equalsIgnoreCase("sift")) {
                    if (dbNSFP == null) dbNSFP = new dbNSFP_Features();
                    dbNSFP.getSIFT(v);
                }
            }
        }

        return true; // if you got here then you had no problems in this function
    }


    // Function applies the user's filter to this variant. Calling this function assumes VariantInfo::fetchFeature()
    // ran successfully. The curTS_ID is only needed if the user requested to evaluate the 'EFF' info field.
    public boolean filter(String jexl_filter_str, String curTS_ID) {

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

        JexlEngine jexl = new JexlBuilder().create(); // Create a jexl engine
        JexlExpression expr = jexl.createExpression(jexl_filter_str); // define the expression you want to test/use

        // Create a MapContext object and populate it with the variables that are in VariantInfo objects
        JexlContext jc = new MapContext();

        if(dbNSFP != null) prepContextMap(dbNSFP, jc, curTS_ID);
        if(ESP != null) prepContextMap(ESP, jc, curTS_ID);
        if(EFF != null) prepContextMap(EFF, jc, curTS_ID);

        boolean retVal = false;
        retVal = (Boolean) expr.evaluate(jc);

        return retVal;
    }



    // Function to populate the jexlContext map with the variables the user has elected to filter on
    private void prepContextMap(Object o, JexlContext jc, String curTS_ID) {


        String dataType = o.getClass().getSimpleName();

        if(dataType.equalsIgnoreCase("dbNSFP_Features")) {
            jc.set("GERP", dbNSFP.GERP__RS);

            if(dbNSFP.mutationTaster.length() > 0) jc.set("MUTATIONTASTER", dbNSFP.mutationTaster);
            if(dbNSFP.polyphen2_hvar.length() > 0) jc.set("POLYPHEN2_VAR", dbNSFP.polyphen2_hvar);
            if(dbNSFP.SIFT.length() > 0) jc.set("SIFT", dbNSFP.SIFT);
        }

        if(dataType.equalsIgnoreCase("ESP_Features")) {
            jc.set("ESP_MAX_AA_EA", ESP.ESP_MAX_AA_AE);
            jc.set("ESP_MAF", ESP.ESP_MAF);
        }

        if(dataType.equalsIgnoreCase("EFF_Features")) {
            String jj = EFF.findTS(curTS_ID);
            jc.set("EFF", EFF.findTS(curTS_ID));
        }
    }

}
