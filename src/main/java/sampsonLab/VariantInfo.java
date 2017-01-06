package sampsonLab;

import htsjdk.variant.variantcontext.VariantContext;

import java.util.Map;


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

            if(feat.equalsIgnoreCase("IS_LOF")) {
                EFF = new EFF_Features( v );
            }
            if(feat.equalsIgnoreCase("MutationTaster")) {
                if(dbNSFP == null) dbNSFP = new dbNSFP_Features();
                dbNSFP.mutationTaster = (String) v;
            }

            if(feat.equalsIgnoreCase("GERP++")) {
                if(dbNSFP == null) dbNSFP = new dbNSFP_Features();
                dbNSFP.GERP__RS = Double.valueOf( (String) v );
            }

            if(feat.equalsIgnoreCase("Polyphen2_Hvar")) {
                if(dbNSFP == null) dbNSFP = new dbNSFP_Features();
                dbNSFP.polyphen2_hvar = (String) v;
            }

            if(feat.equalsIgnoreCase("sift")) {
                if(dbNSFP == null) dbNSFP = new dbNSFP_Features();
                dbNSFP.getSIFT(v);
            }
        }

        return true; // if you got here then you had no problems in this function
    }


    // Function applies the user's filter to this variant
    public boolean filter(String jexl_filter_str) {


        return true;
    }
}
