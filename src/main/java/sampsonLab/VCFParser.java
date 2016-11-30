package sampsonLab;

import htsjdk.variant.vcf.VCFFileReader;

import java.io.File;
import java.util.HashMap;


/**
 * Created by dfermin on 11/30/16.
 */
public class VCFParser {

    static VCFFileReader vcfr = null;
    static File inputVCF = null, inputVCF_tabix = null;

    VCFParser(File vcf, File tbi) {
        inputVCF = vcf;
        inputVCF_tabix = tbi;
    }

    public void parse(HashMap<String, Gene> geneMap) {
        // This command only works if you have the tabix tbi file for the input VCF
        VCFFileReader vcfr = new VCFFileReader(inputVCF, inputVCF_tabix);

        // Get the genomic coordinates for the genes of interest.


    }

}
