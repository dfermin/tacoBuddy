package sampsonLab;

import com.google.common.base.Joiner;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.jexl3.*;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.TreeSet;

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
    public double sample_MAF; // minor allele frequency for the samples

    public HashMap<String, Genotype_Feature> genotypeMap; // k = sampleID, v = genotype_feature object for this subject
    public TreeSet<String> candPatients;

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
        this.sample_MAF = -1.0;

        // construct snp_id string
        snp_id = chr + ":" + String.valueOf(pos) + ":" + REF + ">" + ALT;

        ESP = new ESP_Features();
        dbNSFP = new dbNSFP_Features();

        genotypeMap = new HashMap<String, Genotype_Feature>();
        candPatients = new TreeSet<String>();
    }


    public String getID() { return snp_id; }
    public void setSvmProb(double d) { this.svmProb = d; }


    /*********************************************************************************************/
    // Function records the requested data from the attributes of the variant object
    public boolean fetchFeature(String feat, VariantContext vc) {

        if(feat.startsWith("ESP_")) {
            ESP = new ESP_Features();

            String aa = vc.getAttributeAsString("ESP_AA_AC", ".").replaceAll("\\[", "").replaceAll("\\]", "");
            String ea = vc.getAttributeAsString("ESP_EA_AC", ".").replaceAll("\\[", "").replaceAll("\\]", "");
            String tac = vc.getAttributeAsString("ESP_TAC", ".").replaceAll("\\[", "").replaceAll("\\]", "");
            ESP.ESP_AA_AC = ESP.calcMAF(aa);
            ESP.ESP_EA_AC = ESP.calcMAF(ea);
            ESP.ESP_MAX_AA_EA = Math.max(ESP.ESP_AA_AC, ESP.ESP_EA_AC);
            ESP.ESP_MAF = ESP.calcMAF(tac);
        }

        // Search for this feature in the VariantInfoNameMap
        // If you get a match, record the relevant data for that feature
        if( tacoBuddy.VCF_Info_Name_Map.inMap(feat.toUpperCase()) ) {
            String k = tacoBuddy.VCF_Info_Name_Map.getValue(feat.toUpperCase());

            String v_str = vc.getAttributeAsString(k, "#NULL");
            v_str = v_str.replaceAll("\\[", "").replaceAll("\\]", "");

            if( !v_str.equalsIgnoreCase("#NULL") ) {

                if (feat.equalsIgnoreCase("EFF")) {
                    EFF = new EFF_Features(v_str);
                }
                if (feat.equalsIgnoreCase("MutationTaster")) {
                    if (dbNSFP == null) dbNSFP = new dbNSFP_Features();
                    dbNSFP.mutationTaster = formatStringForJexl(v_str);
                }

                if (feat.equalsIgnoreCase("GERP")) {
                    if (dbNSFP == null) dbNSFP = new dbNSFP_Features();
                    dbNSFP.GERP__RS = Double.valueOf(v_str);
                }

                if (feat.equalsIgnoreCase("Polyphen2_Hvar")) {
                    if (dbNSFP == null) dbNSFP = new dbNSFP_Features();
                    dbNSFP.polyphen2_hvar = formatStringForJexl(v_str);
                }

                if (feat.equalsIgnoreCase("sift")) {
                    if (dbNSFP == null) dbNSFP = new dbNSFP_Features();
                    dbNSFP.SIFT = formatStringForJexl(v_str);
                }
            }
        }

        return true; // if you got here then you had no problems in this function
    }


    /*********************************************************************************************/
    // Function applies the user's filter to this variant. Calling this function assumes VariantInfo::fetchFeature()
    // ran successfully. The curTS_ID is only needed if the user requested to evaluate the 'EFF' info field.
    public boolean filter(String jexl_filter_str) {

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

        JexlEngine jexl = new JexlBuilder().cache(512).strict(true).silent(false).create(); // Create a jexl engine
        JexlExpression expr = jexl.createExpression(jexl_filter_str); // define the expression you want to test/use

        // Create a MapContext object and populate it with the variables that are in VariantInfo objects
        JexlContext jc = new MapContext();

        if(dbNSFP != null) prepContextMap(dbNSFP, jc);
        if(ESP != null) prepContextMap(ESP, jc);
        if(EFF != null) prepContextMap(EFF, jc);

        boolean retVal = false;
        retVal = (Boolean) expr.evaluate(jc);

        if(retVal) passedFilter = 1;

        return retVal;
    }


    /*********************************************************************************************/
    // Function to populate the jexlContext map with the variables the user has elected to filter on
    private void prepContextMap(Object o, JexlContext jc) {

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
            jc.set("EFF", EFF.returnJexlArray());
        }
    }

    /*********************************************************************************************/
    // Function formats the given string so that Jexl will view it as an array.
    // This function assumes the data passed to it is comma separated
    private String formatStringForJexl(String inputStr) {
        String ret = "";
        ArrayList<String> A = new ArrayList<String>();
        for(String s : inputStr.split(",")) {
            String tmp = "'" + s.trim() + "'";
            A.add(tmp);
        }
        ret = "[" + Joiner.on(",").join(A) + "]";
        return ret;
    }

    /*********************************************************************************************/
    // Function returns the requested features for this variant.
    public String returnSummaryString(String geneId, ArrayList<String> PL, String filterType) {

        int PRECISION = 4;
        String ret = "";
        HashMap<String, String> returnValueMap = new HashMap<String, String>();

        for(String feat : globalFunctions.featureSet) {

            returnValueMap.put(feat, "#NULL"); // Initialize map and update it as data is encountered

            // Special case where you have to iterate over the transcripts.
            if (feat.equalsIgnoreCase("EFF")) {
                // output is too long so skip this case
                continue;
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
                returnValueMap.put(feat, dbl2str(dbNSFP.GERP__RS, PRECISION));
                continue;
            }

            if(feat.equalsIgnoreCase("ESP_AA_AC")) {
                returnValueMap.put(feat, dbl2str(ESP.ESP_AA_AC, PRECISION) );
                continue;
            }
            if(feat.equalsIgnoreCase("ESP_EA_AC")) {
                returnValueMap.put(feat, dbl2str(ESP.ESP_EA_AC, PRECISION) );
                continue;
            }
            if(feat.equalsIgnoreCase("ESP_MAF")) {
                returnValueMap.put(feat, dbl2str(ESP.ESP_MAF, PRECISION) );
                continue;
            }
            if(feat.equalsIgnoreCase("ESP_MAX_AA_EA")) {
                returnValueMap.put(feat, dbl2str(ESP.ESP_MAX_AA_EA, PRECISION) );
                continue;
            }
        }

        // Now construct the return string based upon the data in returnValueMap.
        ArrayList<String> dataList = new ArrayList<String>();
        for(String feat : globalFunctions.featureSet) {
            String tmp = returnValueMap.get(feat);
            if(tmp.matches("\\['.'\\]")) tmp = tmp.replaceAll("\\['", "").replaceAll("'\\]", "");
            tmp = tmp.replaceAll("'", "");
            dataList.add(tmp);
        }

        // Append the Sample minor allele frequency (sample_MAF) to the dataList
        dataList.add(dbl2str(sample_MAF, 2));

        String data_line = Joiner.on("\t").join(dataList) + "\n";


        TreeSet<String> patientLines_toMerge = new TreeSet<String>();
        for(String sampleId : PL) {
            Genotype_Feature gf = this.genotypeMap.get(sampleId);

            // We generally aren't interested in people who are homologous for the reference alleles in domiant genes so skip them
            if(gf.genotype_word.equalsIgnoreCase("HOM") && filterType.equalsIgnoreCase("DOM")) continue;

            String tmp = gf.getInfo() + "\t" + filterType + "\t" + geneId + "\t" + snp_id + "\t" + data_line;
            patientLines_toMerge.add(tmp);
        }

        ret = Joiner.on("").join(patientLines_toMerge);

        return ret;
    }

    /*********************************************************************************************/
    // Function adds the observed genotype information for all the subject for the current variant
    public void add(VariantContext vc) {

        ArrayList<String> sampleNames = (ArrayList<String>) vc.getSampleNamesOrderedByName();
        for(String s : sampleNames) {
            Genotype G = vc.getGenotype(s);
            Genotype_Feature gf = new Genotype_Feature(s, G);
            genotypeMap.put(s, gf);
        }
    }


    /*********************************************************************************************/
    // Function returns true if at least 1 subject in the genotypeMap meets the filtering requirements
    public boolean hasCandidateSubjects(String filterType) {
        boolean ret = false;

        for(Genotype_Feature g : genotypeMap.values()) {

            // For dominant genes (DOM) you need to be heterozygous or homologous-domiant
            if(filterType.equalsIgnoreCase("DOM")) {
                if( g.genotype_word.equalsIgnoreCase("HET") ||
                     g.genotype_word.equalsIgnoreCase("HOM") ) candPatients.add(g.sampleID);
            }

            // For recessive genes (REC) you need to be homologous-recessive
            if(filterType.equalsIgnoreCase("REC")) {
                if( g.genotype_word.equalsIgnoreCase("HOM_ALT") ) candPatients.add(g.sampleID);
            }
        }

        if(candPatients.size() > 0) ret = true;

        return ret;
    }

    /*********************************************************************************************/
    // Function returns the minor allele frequency (MAF) for all the subjects in the VCF file.
    // The formula is: 1 - sum_over_i(genotype[0,1,2]) / (2*N)
    // Where i = a subject, and N = number of subjects
    public void calcSampleAF() {

        double sum = 0;
        double n = 0;
        for(Genotype_Feature g : genotypeMap.values()) {
            if(g.genotype_DNA_str.equalsIgnoreCase(".")) continue; // the variant wasn't called for this person so skip this instance
            n++;
            sum += g.genotypeInt;
        }

        double N = 2.0 * n;
        sample_MAF = (1.0 - (sum / N)) * 100; // generally reported as a percent
    }


    /*********************************************************************************************/
    // Function returns the information for the patients in 'candPatients'
    public ArrayList<String> getPatientData() {
        ArrayList<String> ret = new ArrayList<String>(this.candPatients);
        return ret;
    }


    /*********************************************************************************************/
    // Manually round the given double to a specific number of digits and return a String with it's value
    private static String dbl2str(double d, int precision) {
        String t = "%." + Integer.toString(precision) + "f";

        String compStr = "0.";
        for(int i = 0; i < precision; i++) compStr += "0";

        // If after truncation 't' is zero, report 0 instead
        String ret = String.format(t,d);
        if(ret.equalsIgnoreCase(compStr)) ret = "0";

        return ret;
    }
}
