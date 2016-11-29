/**
 * Created by dfermin on 11/29/16.
 */

package sampsonLab;

import java.io.File;
import java.io.IOException;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.*;



public class tacoBuddy {

    static File inputVCF;
    static File inputVCF_tabix;

    public static void main(String[] args) throws IOException {

        if(args.length < 1) {
            System.err.println("\nUSAGE: java -jar tacoBuddy.jar <src_vcf>\n");
            return;
        }

        // Open the VCF file given by the user at the command line.
        // This code assumes you have have a tabix *.tbi file to go with the input VCF file.
        inputVCF = new File(args[0]);
        inputVCF_tabix = new File( new String(args[0] + ".tbi") );

        if( !inputVCF_tabix.exists() ) {
            System.err.println(
                    "\nERROR! Unable to find " + inputVCF_tabix.getName() + " needed. " +
                    "If you don't have it first create it with:\n" +
                    "tabix " + inputVCF.getName() + "\n"
            );
        }

        // This command only works if you have the tabix tbi file for the input VCF
        VCFFileReader vcfr = new VCFFileReader(inputVCF, inputVCF_tabix);


        CloseableIterator<VariantContext> it = vcfr.iterator();
        while(it.hasNext()) {
            VariantContext vc = it.next();
            String contig = vc.getContig(); // this should be the chromosome
            int pos = vc.getEnd(); // position of the variant in the genome
            String ref = vc.getReference().getDisplayString(); // get the reference Allele nucleotide
            String alt = vc.getAltAlleleWithHighestAlleleCount().getDisplayString();

//            List<Allele> AL = vc.getAlleles();
//            Iterator<Allele> Ait = AL.iterator();
//            Allele A = Ait.next();
//            String Astr = A.getDisplayString();

            System.out.print("-------------------------------------------\n");
            System.out.print("chr" + contig + ":" + pos + ref + ">" + alt + "\t" + vc.getID() + "\n");

            /*
            Map<String, Object> M = vc.getAttributes();
            for (String k: M.keySet()) {
                 System.out.print(k + "\n");
            }
            */
            System.out.print("-------------------------------------------\n");
            break;
        }

        vcfr.close();
    }

}
