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
    public int GQ;
    public int totReadDepth, refReadDepth, altReadDepth;
    public String genotype_word; // HOM, HET or HOM_ALT

    public Genotype_Feature(String id, Genotype G) {
        this.sampleID = id;
        this.genotype_DNA_str = "";
        this.genotype_INT_str = "";
        this.genotypeInt = -1;
        this.GQ = -1;
        this.genotype_word = "#NULL";
        this.totReadDepth = 0;
        this.refReadDepth = 0;
        this.altReadDepth = 0;

        if(G.isCalled()) {
            this.genotype_DNA_str = G.getGenotypeString();

            if (G.isHet()) {
                this.genotype_word = "HET";
                this.genotypeInt = 1;
            }
            else if (G.isHomRef()) { // homozygous reference
                this.genotype_word = "HOM";
                this.genotypeInt = 0;
            }
            else if (G.isHomVar()) { // homozygous alternative
                this.genotype_word = "HOM_ALT";
                this.genotypeInt = 2;
            }

            if (G.hasGQ()) this.GQ = G.getGQ();

            if (G.hasAD()) {
                this.refReadDepth = G.getAD()[0];
                this.altReadDepth = G.getAD()[1];
                this.totReadDepth = this.refReadDepth + this.altReadDepth;
            }

            for (Allele a : G.getAlleles()) {
                this.genotype_INT_str += (a.isReference() ? "1" : "0") + "|";
            }
            this.genotype_INT_str = this.genotype_INT_str.substring(0, this.genotype_INT_str.length() - 1);
        }
        else {
            this.genotype_DNA_str = ".";
            this.genotype_INT_str = ".";
        }
    }

    /*----------------------------------------------------------------------------------------------
    ** Function returns the read depth information for this genotype feature
     */
    public String getReadCount() {
        String ret = "";
        ret = "T=" + Integer.toString(this.totReadDepth) + "/" +
              "r=" + Integer.toString(this.refReadDepth) + "/" +
              "a=" + Integer.toString(this.altReadDepth);
        return(ret);
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
