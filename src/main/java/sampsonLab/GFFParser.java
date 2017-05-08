package sampsonLab;

import com.google.common.collect.HashMultimap;

import java.io.*;
import java.util.HashSet;
import java.util.zip.GZIPInputStream;

/**
 * Created by dfermin on 5/3/17.
 */


/*---------------------------------------------------------------------------------------------------------------------
** This class is designed to parse GENCODE formatted GFF3 files. It specifically assumes that's the kind of GFF3 file
** you are going to feed it.
*--------------------------------------------------------------------------------------------------------------------*/
public class GFFParser {

    public File srcFile;


    public GFFParser(File F) {
        srcFile = F;
    }

    public void parse(globalFunctions GF) throws IOException {

        if(GF.genesREC.size() > 0) GF.REC_geneMap = HashMultimap.create(); // k = geneName, v = list of transcripts
        if(GF.genesDOM.size() > 0) GF.DOM_geneMap = HashMultimap.create(); // k = geneName, v = list of transcripts

        if( (GF.genesDOM.size() == 0) && (GF.genesREC.size() == 0)) GF.ALL_geneMap = HashMultimap.create();

        BufferedReader br = null;

        System.err.print("\nParsing " + srcFile.getName() + "\n");

        if(srcFile.getName().endsWith(".gz")) {
            GZIPInputStream gzip = new GZIPInputStream(new FileInputStream(srcFile.getAbsoluteFile()));
            br = new BufferedReader(new InputStreamReader(gzip));
        }
        else {
            FileReader fr = new FileReader(srcFile);
            br = new BufferedReader(fr);
        }


        Transcript curTranscript = null;
        int geneStart = -1, geneEnd = -1;
        char strand = '.';
        String geneID = null, geneType = null, geneName = null, chr = null, line = null;
        boolean isGencodeGFF = false;

        while((line = br.readLine()) != null) {
            if(line.equalsIgnoreCase("#provider: GENCODE")) {
                isGencodeGFF = true;
                continue;
            };

            if(line.startsWith("#")) continue;

            if(!isGencodeGFF) {
                System.err.println("\nERROR: " + srcFile.getName() + " is not a Gencode-generated GFF3 file.\n" +
                "Please download a valid Gencode-gff3 file from https://www.gencodegenes.org/releases/current.html\n");
                System.exit(0);
            }


            String[] data = line.split("\t");
            if(data[2].equalsIgnoreCase("gene")) {

                chr = data[0];
                geneStart = Integer.parseInt(data[3]);
                geneEnd = Integer.parseInt(data[4]);
                strand = data[6].charAt(0);


                String info[] = data[8].split(";");
                for (String s : info) {
                    if (s.startsWith("ID=")) geneID = s.substring(3).replaceAll("\\..+$", "");  // remove version number from ENSG id
                    if (s.startsWith("gene_type=")) geneType = s.substring(10);
                    if (s.startsWith("gene_name=")) geneName = s.substring(10);
                }
            }

            if(data[2].equalsIgnoreCase("transcript")) {

                // Check to see if curTranscript is null, if it isn't you need to store this variable
                // before you can continue;
                if(curTranscript != null) {
                    curTranscript.calcCDSlength();
                    if(GF.genesDOM.contains(curTranscript.getGeneName())) GF.DOM_geneMap.put(curTranscript.getGeneName(), curTranscript);
                    else if(GF.genesREC.contains(curTranscript.getGeneName())) GF.REC_geneMap.put(curTranscript.getGeneName(), curTranscript);
                    else GF.ALL_geneMap.put(curTranscript.getGeneName(), curTranscript);
                }
                curTranscript = null;

                int start = Integer.parseInt(data[3]);
                int end = Integer.parseInt(data[4]);

                String localGeneName = null, transID = null;

                String info[] = data[8].split(";");
                for (String s : info) {
                    if (s.startsWith("ID=")) transID = s.substring(3).replaceAll("\\..+$", ""); // remove version number from ENST id
                    if (s.startsWith("gene_name=")) localGeneName = s.substring(10);
                }

                if (!localGeneName.equalsIgnoreCase(geneName)) {
                    System.err.print("\nERROR! (transcript): " + localGeneName + " != " + geneName + " for current line of srcGFF3:\n" +
                            line + "\n\n");
                    System.exit(0);
                }

                // Create a new transcript entry
                curTranscript = new Transcript(chr,
                        geneStart, geneEnd,
                        strand, geneID, geneType, localGeneName,
                        start, end, transID);
            }


            if(data[2].equalsIgnoreCase("exon")) {
                int start = Integer.parseInt(data[3]);
                int end = Integer.parseInt(data[4]);

                String transID = null, exonID = null;
                int exonNum = 0;
                String info[] = data[8].split(";");
                for (String s : info) {
                    if (s.startsWith("Parent=")) transID = s.substring(7).replaceAll("\\..+$", "");
                    if (s.startsWith("exon_id=")) exonID = s.substring(8).replaceAll("\\..+$", "");
                    if (s.startsWith("exon_number=")) exonNum = Integer.parseInt(s.substring(12));
                }

                if(!transID.equalsIgnoreCase(curTranscript.getTranscriptID())) {
                    System.err.print("\nERROR! (exon): " + transID + " != " + curTranscript.getTranscriptID() + " for current line of srcGFF3:" +
                            line + "\n\n");
                    System.exit(0);
                }

                Exon E = new Exon(exonID,
                        curTranscript.getGeneName(),
                        curTranscript.getChrom(),
                        start, end, exonNum
                );

                curTranscript.addExon(E);
                E = null;
            }
        }
        br.close();

        if(GF.genesDOM.size() > 0) System.err.println("DOM gene map size: " + GF.DOM_geneMap.asMap().size());
        if(GF.genesREC.size() > 0) System.err.println("REC gene map size: " + GF.REC_geneMap.asMap().size());
        if(GF.ALL_geneMap.asMap().size() > 0) System.err.println("Gene map size: " + GF.ALL_geneMap.asMap().size());

    }

}
