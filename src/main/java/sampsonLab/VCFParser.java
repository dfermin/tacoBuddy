package sampsonLab;

import com.google.common.collect.SetMultimap;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

import java.io.File;
import java.util.Set;


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

    public void parse(SetMultimap<String, Transcript> geneMap) {

        // This command only works if you have the tabix tbi file for the input VCF
        VCFFileReader vcfr = new VCFFileReader(inputVCF, inputVCF_tabix);

        // Iterate over the genes in the given geneMap
        for(String geneId: geneMap.keySet()) {

            Set<Transcript> allTS = geneMap.get(geneId);

            for(Transcript curTS : allTS) {

                CloseableIterator<VariantContext> it = vcfr.query(
                        curTS.getTrimmedChrom(),
                        curTS.get_ts_Start(),
                        curTS.get_ts_End());

                while(it.hasNext()) {
                    VariantContext vc = it.next();
                    String chr = vc.getContig();
                    int pos = vc.getEnd();
                    String ref = vc.getReference().getDisplayString(); // get the reference Allele nucleotide
                    String alt = vc.getAltAlleleWithHighestAlleleCount().getDisplayString();
                    String dbSNP_id = vc.getID();

                    System.out.print(
                            curTS.getGeneName() + "\t" + curTS.getTranscriptID() + "\t" +
                            curTS.getChrom() + ":" + curTS.get_ts_Start() + ".." + curTS.get_ts_End() + "\t" +
                            "chr" + chr + ":" + pos + ref + ">" + alt + "\t" + dbSNP_id + "\n"
                    );

                }
            }
        }
    }

}
