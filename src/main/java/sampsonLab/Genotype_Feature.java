package sampsonLab;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;


/**
 * Created by dfermin on 1/16/17.
 */

public class Genotype_Feature extends FeatureClass {
    public String sampleID;
    public String genotype_DNA_str; // genotype represented as nucleotides, example: G/A or T|T
    public String genotype_INT_str; // genotype represented as numbers, example: 0/1, 1/1, 0/0, 1/0
    public int genotypeInt; // 0 = homologous dominant, 1 = heterozygous, 2 = homologous recessive
    public int DP, GQ;
    public String genotype_word; // HOM HET or HOM_ALT

    public Genotype_Feature(String id, Genotype G) {
        this.sampleID = id;
        this.genotype_DNA_str = "";
        this.genotype_INT_str = "";
        this.genotypeInt = -1;
        this.DP = -1;
        this.GQ = -1;
        this.genotype_word = "";

        if(G.isCalled()) {
            this.genotype_DNA_str = G.getGenotypeString();

            if (G.isHet()) this.genotype_word = "HET";
            else if (G.isHomRef()) this.genotype_word = "HOM";
            else if (G.isHomVar()) this.genotype_word = "HOM_ALT";

            if (G.hasAD()) this.DP = G.getDP();
            if (G.hasGQ()) this.GQ = G.getGQ();

            for (Allele a : G.getAlleles()) {
                this.genotype_INT_str += (a.isReference() ? "1" : "0") + "|";
            }
            this.genotype_INT_str = this.genotype_INT_str.substring(0, this.genotype_INT_str.length() - 1);

            if (this.genotype_INT_str.equalsIgnoreCase("0|0")) genotypeInt = 0;
            if ((this.genotype_INT_str.equalsIgnoreCase("1|0")) ||
                    (this.genotype_INT_str.equalsIgnoreCase("0|1"))) genotypeInt = 1;
            if (this.genotype_INT_str.equalsIgnoreCase("1|1")) genotypeInt = 2;
        }
        else {
            this.genotype_DNA_str = ".";
            this.genotype_INT_str = ".";
        }
    }

    /**********************************************************************************************/
    // Function returns tab-delimited string of the information for this person we have collected for
    // this genotype call
    public String getInfo() {
        String ret = "";

        if( !genotype_DNA_str.equalsIgnoreCase(".") )
            ret += this.sampleID + "\t" + this.genotype_word;

        return ret;
    }
}
