/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package proteasixcsdatabase;

import com.sun.media.jai.opimage.AddCRIF;
import java.io.*;
import java.net.MalformedURLException;
import java.net.URL;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.SQLException;
import java.sql.Statement;
import java.text.SimpleDateFormat;
import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import org.w3c.dom.Document;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

/**
 *
 * @author julieklein
 */
public class ProteasixCSdatabase {

    private StringBuilder getHtmlcontent(URL u) throws IOException {
        InputStream is = null;
        DataInputStream dis;
        String s = null;
        StringBuilder htmlcontent = new StringBuilder();

        is = u.openStream();
        dis = new DataInputStream(new BufferedInputStream(is));
        while ((s = dis.readLine()) != null) {
            s = s.replaceAll("\\?", "");
            htmlcontent.append(s);
        }
        is.close();
        return htmlcontent;
    }

    private void getManualentries(LinkedList<CsDatabaseEntry> alldatanotcurated) throws IOException, FileNotFoundException, NumberFormatException {
        //1c. Manual

        BufferedReader bReader = createBufferedreader("//Users/julieklein/Dropbox/ProteasiX/ProteasiX/MC_ManualCSDBVersion20120201.txt");
        String line;
        while ((line = bReader.readLine()) != null) {
            String splitarray[] = line.split("\t");
            String proteaseUni = splitarray[1];
            String proteaseBrenda = splitarray[2];
            String substrateUni = splitarray[5];
            String p1position = splitarray[8];
            int ip1position = Integer.parseInt(p1position);
            String p1primeposition = splitarray[9];
            int ip1primeposition = Integer.parseInt(p1primeposition);
            String p1sequence = splitarray[10];
            String p1primesequence = splitarray[11];
            String externallink = splitarray[12];
            String pmid = splitarray[13];

            CsDatabaseEntry csdatabase = new CsDatabaseEntry();
            ProteaseDatabaseEntry proteasedatabase = new ProteaseDatabaseEntry();
            SubstrateDatabaseEntry substratedatabase = new SubstrateDatabaseEntry();


            SimpleDateFormat format = new SimpleDateFormat("dd-MM-yyyy");
            Calendar originalDate = Calendar.getInstance();
            String dateString = format.format(originalDate.getTime());
            System.out.println(dateString);
            csdatabase.setCreation_Date(dateString);

            csdatabase.setPMID(pmid);
            csdatabase.setExternal_Link(externallink);

            csdatabase.setP1_Position(ip1position);
            csdatabase.setP1prime_Position(ip1primeposition);
            csdatabase.setP1_Sequence(p1sequence);
            csdatabase.setP1prime_Sequence(p1primesequence);

            substratedatabase.setS_UniprotID(substrateUni);
            String UniprotUrlsub = "http://www.uniprot.org/uniprot/" + substrateUni + ".xml";
            NodeList entriesmanualsub = getEntries("/uniprot/entry", parseUniprot(UniprotUrlsub));

            for (int i = 0; i < entriesmanualsub.getLength(); i++) {
                getUniSubstratepproteinname(entriesmanualsub, i, substratedatabase);
                getUniSubstrategenename(entriesmanualsub, i, substratedatabase);
                String substrateTaxon = getUniSubstratetaxonomy(entriesmanualsub, i, substratedatabase);
                proteasedatabase.setP_Taxon(substrateTaxon);
                String sequence = getUniSubstratesequence(entriesmanualsub, i, substratedatabase);

                int sequencelength = sequence.length();

                if (ip1position - 4 > 0 || ip1position - 4 == 0) {
                    if (ip1primeposition + 3 == sequencelength || ip1primeposition + 3 < sequencelength) {
                        String cleavagesite = sequence.substring(ip1position - 4, ip1primeposition + 3);
                        csdatabase.setCleavagesiteseaquence(cleavagesite);
                    } else {
                        String cleavagesite = sequence.substring(ip1position - 4, sequencelength);
                        csdatabase.setCleavagesiteseaquence(cleavagesite);
                    }
                } else {
                    if (ip1primeposition + 3 == sequencelength || ip1primeposition + 3 < sequencelength) {
                        String cleavagesite = sequence.substring(0, ip1primeposition + 3);
                        csdatabase.setCleavagesiteseaquence(cleavagesite);
                    } else {
                        String cleavagesite = sequence.substring(0, sequencelength);
                        csdatabase.setCleavagesiteseaquence(cleavagesite);
                    }
                }

            }


            proteasedatabase.setP_UniprotID(proteaseUni);
            proteasedatabase.setP_EC_Number(proteaseBrenda);
            String UniprotUrlprot = "http://www.uniprot.org/uniprot/" + proteaseUni + ".xml";
            NodeList entriesmanualprot = getEntries("/uniprot/entry", parseUniprot(UniprotUrlprot));
            for (int i = 0; i < entriesmanualprot.getLength(); i++) {
                getUniProteasegenename(entriesmanualprot, i, proteasedatabase);
                getUniProteasepproteinname(entriesmanualprot, i, proteasedatabase);
            }
            csdatabase.setComment("-");
            csdatabase.setCuration_Status("-");
            csdatabase.setProtease(proteasedatabase);
            csdatabase.setSubstrate(substratedatabase);
            alldatanotcurated.add(csdatabase);
        }
    }

    private Matcher getPatternmatcher(String expression, String string) {
        Pattern p = Pattern.compile(expression,
                Pattern.DOTALL | Pattern.UNIX_LINES | Pattern.MULTILINE);
        Matcher matcher = p.matcher(string);
        return matcher;
    }

    private String mapProteasetoLibrairy(String commentS, String proteaseTaxon, String proteaseName, CsDatabaseEntry csdatabase, ProteaseDatabaseEntry proteasedatabase) throws IOException {
        String commentP;
        commentP = commentS + "; Check Protease Symbol and Accession; add to Substrate Librairy";
        if (proteaseTaxon.equalsIgnoreCase("Homo Sapiens")) {
            BufferedReader bReader = createBufferedreader("/Users/julieklein/Dropbox/ProteasiX/LIBRAIRIES/ProteaseHSALibrairy.txt");
            commentP = getProteaseInformation(bReader, proteaseName, csdatabase, proteasedatabase, commentS);
        } else {
            if (proteaseTaxon.equalsIgnoreCase("Rattus Norvegicus")) {
                BufferedReader bReader = createBufferedreader("/Users/julieklein/Dropbox/ProteasiX/LIBRAIRIES/ProteaseRNOLibrairy.txt");
                commentP = getProteaseInformation(bReader, proteaseName, csdatabase, proteasedatabase, commentS);
            } else {
                if (proteaseTaxon.equalsIgnoreCase("Mus Musculus")) {
                    BufferedReader bReader = createBufferedreader("/Users/julieklein/Dropbox/ProteasiX/LIBRAIRIES/ProteaseMMULibrairy.txt");
                    commentP = getProteaseInformation(bReader, proteaseName, csdatabase, proteasedatabase, commentS);
                }
            }
        }
        return commentP;
    }

    private String putSplittedhtmlintostringbuilder(Matcher splithtml) {
        StringBuilder sbd = new StringBuilder();
        while (splithtml.find()) {
            String entry = splithtml.group(1);
            entry = entry + "\n" + "******************************" + "\n";
            //System.out.println(entry);
            sbd.append(entry);
        }
        String splittedentry = sbd.toString();
        return splittedentry;
    }

    private String getProteaseTaxon(Matcher patternProteaseTaxon) {
        String proteaseTaxon = null;
        if (patternProteaseTaxon.find()) {
            proteaseTaxon = patternProteaseTaxon.group(2);
            proteaseTaxon = proteaseTaxon.trim();
        } else {
            proteaseTaxon = "n.d.";
        }
        return proteaseTaxon;
    }

    private String getProteaseNameSymbolId(Matcher patternProteaseName, CsDatabaseEntry csdatabase, ProteaseDatabaseEntry proteasedatabase, String proteaseTaxon, String commentS) throws IOException {
        String commentP = null;
        if (patternProteaseName.find()) {
            String proteaseName = patternProteaseName.group(1);
            proteaseName = proteaseName.trim();
            proteaseName = proteaseName.replaceAll(",", "");
            proteaseName = proteaseName.replaceAll(";", "");
            commentP = mapProteasetoLibrairy(commentS, proteaseTaxon, proteaseName, csdatabase, proteasedatabase);

        } else {
            String proteaseName = "n.d.";
            String proteaseSymbol = "n.d.";
            String proteaseUniprot = "n.d";
            String proteaseBrenda = "n.d.";
            commentP = commentS + ";-";
            proteasedatabase.setP_NL_Name(proteaseName);
            proteasedatabase.setP_Symbol(proteaseSymbol);
            proteasedatabase.setP_UniprotID(proteaseUniprot);
            proteasedatabase.setP_EC_Number(proteaseBrenda);
            System.out.println(proteaseName);
            System.out.println(proteaseSymbol);
            System.out.println(proteaseUniprot);
            System.out.println(proteaseBrenda);
        }
        return commentP;
    }

    private BufferedReader createBufferedreader(String datafilename) throws FileNotFoundException {
        BufferedReader bReader = new BufferedReader(
                new FileReader(datafilename));
        return bReader;

    }

    private String getProteaseInformation(BufferedReader bReader, String proteaseName, CsDatabaseEntry csdatabase, ProteaseDatabaseEntry proteasedatabase, String commentS) throws IOException {
        String line;
        String commentP = null;
        while ((line = bReader.readLine()) != null) {
            String splitarray[] = line.split("\t");
            String naturallanguage = splitarray[0];
            naturallanguage = naturallanguage.replaceAll("\"", "");
            if (naturallanguage.equals(proteaseName)) {
                String proteaseSymbol = splitarray[1];
                proteaseSymbol = proteaseSymbol.replaceAll("sept-0", "SEPT");
                String proteaseUniprot = splitarray[2];
                String proteaseBrenda = splitarray[3];
                String UniprotURL = "http://www.uniprot.org/uniprot/" + proteaseUniprot + ".xml";
                NodeList entries = getEntries("/uniprot/entry", parseUniprot(UniprotURL));
                for (int i = 0; i < entries.getLength(); i++) {
                    getUniProteasepproteinname(entries, i, proteasedatabase);
                    String genename = getUniProteasegenename(entries, i, proteasedatabase);
                }
                commentP = commentS + ";-";
                proteasedatabase.setP_UniprotID(proteaseUniprot);
                proteasedatabase.setP_EC_Number(proteaseBrenda);
                csdatabase.setProtease(proteasedatabase);
                System.out.println(proteaseUniprot);
                System.out.println(proteaseBrenda);

            }

        }
        return commentP;
    }

    private String getSubstrateNameSymbolId(Matcher patternSubstrateName, SubstrateDatabaseEntry substratedatabase, String entry) throws IOException {
        String commentS = null;
        if (patternSubstrateName.find()) {
            String Substratename = patternSubstrateName.group(1);
            Substratename = Substratename.trim();
            Substratename = Substratename.replaceAll(",", "");
            Substratename = Substratename.replaceAll(";", "");
            commentS = "Check Substrate Symbol and Accession; add to Substrate Librairy";
            Matcher patternSubstrateSymbol = getPatternmatcher("Substrate[^<]+</th>[^<]+<td>[^<]+<table>[^<]+<tr>[^<]+<th\\s+class=\"th3\">Definition:</th>[^<]+<td><b>"
                    + "[^<]+"
                    + "</b></td>[^<]+<tr>[^<]+<th\\s+class=\"th3\">Symbol:</th>"
                    + "([^<]+<td><b>)?([^<]+)", entry);
            getSubstrateSymbol(patternSubstrateSymbol, substratedatabase);
            Matcher patternSubstrateAccession = getPatternmatcher("UniProt\\s+Accession:</th>[^<]+<td><a\\s+href\\s+=\\s+\"[^\"]+\"\\s+target=\"[^\"]+\">"
                    + "([^<]+)", entry);
            getSubstrateAccession(patternSubstrateAccession, substratedatabase);

            BufferedReader bReader = createBufferedreader("/Users/julieklein/Dropbox/ProteasiX/LIBRAIRIES/SubstrateHSALibrairy.txt");
            String line;
            while ((line = bReader.readLine()) != null) {
                String splitarray[] = line.split("\t");
                String naturallanguage = splitarray[1];
                naturallanguage = naturallanguage.replaceAll("\"", "");
                if (naturallanguage.equals(Substratename)) {
                    String Substratesymbol = splitarray[0];
                    Substratesymbol = Substratesymbol.replaceAll("sept-0", "SEPT");
                    String Substrateaccession = splitarray[2];
                    String UniprotURL = "http://www.uniprot.org/uniprot/" + Substrateaccession + ".xml";
                    NodeList entries = getEntries("/uniprot/entry", parseUniprot(UniprotURL));
                    for (int i = 0; i < entries.getLength(); i++) {
                        getUniSubstratepproteinname(entries, i, substratedatabase);
                        String genename = getUniSubstrategenename(entries, i, substratedatabase);
                    }
                    System.out.println(Substrateaccession);

                    substratedatabase.setS_UniprotID(Substrateaccession);
                    commentS = "-";
                    System.out.println(commentS);

                }
            }

        } else {
            String Substratename = "n.d.";
            String Substratesymbol = "n.d.";
            String Substrateaccession = "n.d";
            substratedatabase.setS_NL_Name(Substratename);
            substratedatabase.setS_Symbol(Substratesymbol);
            substratedatabase.setS_UniprotID(Substrateaccession);
            System.out.println(Substratename);
            System.out.println(Substratesymbol);
            System.out.println(Substrateaccession);
            commentS = "-";
            System.out.println(commentS);
        }
        return commentS;
    }

    private void getSubstrateSymbol(Matcher patternSubstrateSymbol, SubstrateDatabaseEntry substratedatabase) {

        if (patternSubstrateSymbol.find()) {
            String Substratesymbol = patternSubstrateSymbol.group(2);
            substratedatabase.setS_Symbol(Substratesymbol);
            //System.out.println(Substratesymbol);
        } else {
            String Substratesymbol = "n.d.";
            substratedatabase.setS_Symbol(Substratesymbol);
            //System.out.println(Substratesymbol);
        }
    }

    private void getSubstrateAccession(Matcher patternSubstrateAccession, SubstrateDatabaseEntry subtratedatabase) {
        if (patternSubstrateAccession.find()) {
            String accession = patternSubstrateAccession.group(1);
            accession = accession.trim();
            subtratedatabase.setS_UniprotID(accession);
            //System.out.println(accession);
        } else {
            String accession = "n.d.";
            subtratedatabase.setS_UniprotID(accession);
            //System.out.println(accession);
        }
    }

    private void getCleavagesitePosition(Matcher patternCleavagesitePosition, CsDatabaseEntry csdatabase) throws NumberFormatException {
        while (patternCleavagesitePosition.find()) {
            String position = patternCleavagesitePosition.group(1);
            position = position.trim();
            if (position.equalsIgnoreCase("No_information")) {
                int intP1 = 0;
                csdatabase.setP1_Position(intP1);
                int intP1prime = 0;
                csdatabase.setP1prime_Position(intP1prime);
                System.out.println(intP1);
                System.out.println(intP1prime);
            } else {
                //System.out.println(position);
                String positionSplit[] = position.split("-");
                String P1 = positionSplit[0];
//                    String P1prime = positionSplit[1];
                int intP1 = Integer.parseInt(P1);
                csdatabase.setP1_Position(intP1);
//                    int intP1prime = Integer.parseInt(P1prime);
                int intP1prime = intP1 + 1;
                csdatabase.setP1prime_Position(intP1prime);
                System.out.println(intP1);
                System.out.println(intP1prime);
            }
        }
    }

    private String getCleavagesiteSequence(Matcher patternCleavagesiteSequence, CsDatabaseEntry csdatabase) {
        String csSequence = null;
        while (patternCleavagesiteSequence.find()) {
            csSequence = patternCleavagesiteSequence.group(1);
            csSequence = csSequence.trim();
            if (csSequence.equalsIgnoreCase("No_information")) {
                String aaP1 = "?";
                String aaP1prime = "?";
                csdatabase.setP1_Sequence(aaP1);
                csdatabase.setP1prime_Sequence(aaP1prime);
                System.out.println(aaP1);
                System.out.println(aaP1prime);
                csdatabase.setCleavagesiteseaquence("no information");
                System.out.println("no information");

            } else {
                //System.out.println(position);
                String positionSplit[] = csSequence.split("");
                String aaP1 = positionSplit[4];
                String aaP1prime = positionSplit[6];
                csdatabase.setP1_Sequence(aaP1);
                csdatabase.setP1prime_Sequence(aaP1prime);
                System.out.println(aaP1);
                System.out.println(aaP1prime);
                String csSequencenodash = csSequence.replaceAll("-", "");
                csdatabase.setCleavagesiteseaquence(csSequencenodash);
                System.out.println(csSequencenodash);
            }
        }
        return csSequence;
    }

    private void getPmid(Matcher patternPmid, CsDatabaseEntry csdatabase) {
        if (patternPmid.find()) {
            String pmid = patternPmid.group(1);
            pmid = pmid.trim();
            pmid = pmid.replaceAll(",", ";");
            csdatabase.setPMID(pmid);
            System.out.println(pmid);
        } else {
            String pmid = "-";
            csdatabase.setPMID(pmid);
            System.out.println(pmid);
        }
    }

    private void getErrorUnmatched(Matcher patternErrorUnmatched, CsDatabaseEntry csdatabase) {
        if (patternErrorUnmatched.find()) {
            String errormunmatched = patternErrorUnmatched.group(1);
            String error = "Unmatched cleavage site; Cleavage site discarded";
            csdatabase.setCuration_Status(error);
            System.out.println(error);

        } else {
            String error = "-";
            csdatabase.setCuration_Status(error);
            System.out.println(error);
        }

    }

    private Document parseUniprot(String url) {
        ParseUniprot parser = new ParseUniprot();
        Document xml = parser.getXML(url);
        xml.getXmlVersion();
        return xml;
    }

    public static int getLineCount(String url) throws IOException {
        int linecounter = 0;
        BufferedReader urlReader = new BufferedReader(new InputStreamReader(new URL(url).openStream()));
        while (urlReader.readLine() != null) {
            linecounter++;
        }
        System.out.println(linecounter);
        return linecounter;
    }

    private String getUniSubstrategenename(NodeList entries, int i, SubstrateDatabaseEntry substratedatabase) {
        //GET SUBSTRATE GENE NAME using getInformation method
        LinkedList<String> genenamelist = getInformation("./gene/name[@type][1]/text()", entries.item(i));
        String genename = null;
        if (!genenamelist.isEmpty()) {
            genename = genenamelist.getFirst();
            System.out.println(genename);
            substratedatabase.setS_Symbol(genename);
        }
        return genename;
    }

    private void getUniSubstratepproteinname(NodeList entries, int i, SubstrateDatabaseEntry substratedatabase) {
        //GET SUBSTRATE PROTEIN NAME using getInformation method
        LinkedList<String> protnamelist = getInformation("./protein/recommendedName/fullName/text()", entries.item(i));
        String protname = null;
        if (!protnamelist.isEmpty()) {
            protname = protnamelist.getFirst();
            protname = protname.replaceAll(",", "");
            System.out.println(protname);
            substratedatabase.setS_NL_Name(protname);
        }
    }

    private String getUniProteasegenename(NodeList entries, int i, ProteaseDatabaseEntry proteasedatabase) {
        //GET SUBSTRATE GENE NAME using getInformation method
        LinkedList<String> genenamelist = getInformation("./gene/name[@type][1]/text()", entries.item(i));
        String genename = null;
        if (!genenamelist.isEmpty()) {
            genename = genenamelist.getFirst();
            System.out.println(genename);
            proteasedatabase.setP_Symbol(genename);
        }
        return genename;
    }

    private void getUniProteasepproteinname(NodeList entries, int i, ProteaseDatabaseEntry proteasedatabase) {
        //GET SUBSTRATE PROTEIN NAME using getInformation method
        LinkedList<String> protnamelist = getInformation("./protein/recommendedName/fullName/text()", entries.item(i));
        String protname = null;
        if (!protnamelist.isEmpty()) {
            protname = protnamelist.getFirst();
            protname = protname.replaceAll(",", "");
            System.out.println(protname);
            proteasedatabase.setP_NL_Name(protname);
        }
    }

    private String getUniSubstratesequence(NodeList entries, int i, SubstrateDatabaseEntry substratedatabase) {
        //GET PROTSEQUENCE using getInformation method
        String sequence = getInformation("/uniprot/entry/sequence/text()", entries.item(i)).getFirst();
        sequence = sequence.replaceAll("\n", "");
        return sequence;
    }

    private String getUniSubstratetaxonomy(NodeList entries, int i, SubstrateDatabaseEntry substratedatabase) {
        //GET TAXONOMY using getInformation method
        String taxon = getInformation("./organism/dbReference/@id", entries.item(i)).getFirst();
        taxon.trim();
        if (taxon.equalsIgnoreCase("9606")) {
            taxon = "Homo Sapiens";
        }
        System.out.println(taxon);
        substratedatabase.setS_Taxon(taxon);
        return taxon;
    }

    private String getAccession(NodeList entries, int i, SubstrateDatabaseEntry substratedatabase) {
        //GET URL AND GET ACCESSION using getInformation method
        LinkedList<String> idlist = getInformation("./accession[1]/text()", entries.item(i));
        String uniprotid = idlist.getFirst();
        System.out.println(uniprotid);
        substratedatabase.setS_UniprotID(uniprotid);
        return uniprotid;
    }

    private NodeList getEntries(String query, Document xml) {
        XPathUniprot XPather = new XPathUniprot();
        NodeList entrylist = XPather.getNodeListByXPath(query, xml);
        return entrylist;
    }

    private LinkedList<String> getInformation(String query, Node i) {
        XPathNodeUniprot XPathNoder = new XPathNodeUniprot();
        NodeList entrynodelist = XPathNoder.getNodeListByXPath(query, i);
        Loop l1 = new Loop();
        LinkedList<String> information = l1.getStringfromNodelist(entrynodelist);
        return information;
    }

    private NodeList getCsdbentries(String query, Node i) {
        XPathNodeUniprot XPathNoder = new XPathNodeUniprot();
        NodeList csdbentries = XPathNoder.getNodeListByXPath(query, i);
        return csdbentries;
    }

    public ProteasixCSdatabase() throws MalformedURLException, IOException {
        PrintStream csvWriter = null;
        LinkedList<CsDatabaseEntry> alldatanotcurated = new LinkedList<CsDatabaseEntry>();
        LinkedList<CsDatabaseEntry> alldatacurated = new LinkedList<CsDatabaseEntry>();
        java.util.Calendar calendar = java.util.Calendar.getInstance();

        String version = "AVRIL2012";


        //1. Create CsDatabase "; separated file"

        //1a. PMAP

//        File f = new File("//Users/julieklein/Dropbox/ProteasiX/ProteasiX/PMAP/AVRIL2012/");
//        File[] files = f.listFiles();
//        for (File file : files) {
//            String filepath = "file://" + file.getPath() + "/";
//            String htmlcontentmultipleentries = getHtmlcontent(new URL(filepath)).toString();
//
//            Matcher splithtml = getPatternmatcher("(<input\\s+id=\"ballot.*?>Detail</a></td>)", htmlcontentmultipleentries);
//            String htmlsplitted = putSplittedhtmlintostringbuilder(splithtml);
//            Matcher retrievepmapentryid = getPatternmatcher("<td><a\\s+href=\""
//                    + "([^\"]+)"
//                    + "\"[^>]*>", htmlsplitted);
//
//            while (retrievepmapentryid.find()) {
//                String url = retrievepmapentryid.group(1);
//                if (url.equalsIgnoreCase("/relation/show/16398") || url.equalsIgnoreCase("/relation/show/17178") || url.equalsIgnoreCase("/relation/show/17177") || url.equalsIgnoreCase("/relation/show/17074") || url.equalsIgnoreCase("/relation/show/17458") || url.equalsIgnoreCase("/relation/show/17467") || url.equalsIgnoreCase("/relation/show/16083") || url.equalsIgnoreCase("/relation/show/16082") || url.equalsIgnoreCase("/relation/show/16081") || url.equalsIgnoreCase("/relation/show/16080") || url.equalsIgnoreCase("/relation/show/16398") || url.equalsIgnoreCase("/relation/show/16271")) {
//                } else {
//                    url = "http://cutdb.burnham.org" + url;
//                    String entry = getHtmlcontent(new URL(url)).toString();
//                    Matcher patternProteaseTaxon = getPatternmatcher("<div\\s+id=\"protdata\">[^<]+<table>[^<]+<tr>[^<]+<th\\s+class=\"th3\">Definition:</th>"
//                            + "([^<]+<td><b>[^<]+</b></td>)?"
//                            + "[^<]+</tr>[^<]+<tr>[^<]+<th\\s+class=\"th3\">Organism:</th>[^<]+<td><a\\s+href\\s+=\\s+\"[^\"]+\"\\s+target=\"[^\"]+\">"
//                            + "([^<]+)", entry);
//                    String proteaseTaxon = getProteaseTaxon(patternProteaseTaxon);
//
//                    if (!(proteaseTaxon.equalsIgnoreCase("Homo Sapiens")) && !(proteaseTaxon.equalsIgnoreCase("Rattus Norvegicus")) && !(proteaseTaxon.equalsIgnoreCase("Mus Musculus"))) {
//                        continue;
//                    } else {
//
//                        System.out.println("\n" + "******************************" + url);
//
//
//
//                        Matcher patternSubstrateName = getPatternmatcher("Substrate[^<]+</th>[^<]+<td>[^<]+<table>[^<]+<tr>[^<]+<th\\s+class=\"th3\">Definition:</th>[^<]+<td><b>"
//                                + "([^<]+)", entry);
//                        SubstrateDatabaseEntry substratedatabase = new SubstrateDatabaseEntry();
//
//                        String commentS = getSubstrateNameSymbolId(patternSubstrateName, substratedatabase, entry);
//
//                        CsDatabaseEntry csdatabase = new CsDatabaseEntry();
//                        ProteaseDatabaseEntry proteasedatabase = new ProteaseDatabaseEntry();
//                        proteasedatabase.setP_Taxon(proteaseTaxon);
//                        System.out.println(proteaseTaxon);
//                        csdatabase.setSubstrate(substratedatabase);
//                                
//                        Matcher patternProteaseName = getPatternmatcher("<div\\s+id=\"protdata\">[^<]+<table>[^<]+<tr>[^<]+<th\\s+class=\"th3\">Definition:</th>[^<]+<td><b>"
//                                + "([^<]+)", entry);
//                        String commentP = getProteaseNameSymbolId(patternProteaseName, csdatabase, proteasedatabase, proteaseTaxon, commentS);       
//
//                        
//                        csdatabase.setComment(commentP);
//                        System.out.println(commentP);
//                        csdatabase.setExternal_Link(url);
//
//                        Matcher patternCleavagesitePosition = getPatternmatcher("<div\\s+id=\"cleav2\">[^<]+<table>[^<]+<th\\s+class=\"th3\">Position:</th>[^<]+<td>"
//                                + "([^<]+)?", entry);
//                        getCleavagesitePosition(patternCleavagesitePosition, csdatabase);
//
//                        Matcher patternCleavagesiteSequence = getPatternmatcher("<div\\s+id=\"cleav\">[^<]+<table>[^<]+<th\\s+class=\"th3\">Sequence:</td>[^<]+<td>"
//                                + "([^<]+)?", entry);
//                        String csSequence = getCleavagesiteSequence(patternCleavagesiteSequence, csdatabase);
//
//                        Matcher patternPmid = getPatternmatcher("<div\\s+id=\"pubmed\">[^<]+<table>[^<]+<td>[^<]+<a\\s+href=\"[^\"]+\"\\s+target=\"[^\"]+\"\\s+>"
//                                + "([^<]+)", entry);
//                        getPmid(patternPmid, csdatabase);
//
//                        Matcher patternErrorUnmatched = getPatternmatcher("<td><font\\s+color=\"#FF0000\">"
//                                + "(\\*Unmatched)", entry);
//                        getErrorUnmatched(patternErrorUnmatched, csdatabase);
//
//                        SimpleDateFormat format = new SimpleDateFormat("dd-MM-yyyy");
//                        Calendar originalDate = Calendar.getInstance();
//                        String dateString = format.format(originalDate.getTime());
//                        System.out.println(dateString);
//                        csdatabase.setCreation_Date(dateString);
//
//
//                        if (!substratedatabase.getS_UniprotID().equalsIgnoreCase("n.d") && !csdatabase.getP1_Sequence().equalsIgnoreCase("?") && !csdatabase.getP1prime_Sequence().equalsIgnoreCase("?")) {
//                            String motif = csdatabase.getP1_Sequence() + "-" + csdatabase.getP1prime_Sequence();
//                            int motifposition = csSequence.indexOf(motif) + 1;
////                            System.out.println(motif);
////                            System.out.println(motifposition);
//
//                            String csSequencenodash = csSequence.replaceAll("-", "");
//
//                            NodeList entries = getEntries("/uniprot/entry/sequence/text()", parseUniprot("http://www.uniprot.org/uniprot/" + substratedatabase.getS_UniprotID() + ".xml"));
//                            for (int i = 0; i < entries.getLength(); i++) {
//                                String sequence = getUniSubstratesequence(entries, i, substratedatabase);
//
//                                if (sequence.contains(csSequencenodash)) {
//                                    int cleavagesiteposition = sequence.indexOf(csSequencenodash);
//                                    int newp1 = motifposition + cleavagesiteposition;
//                                    int newp1prime = newp1 + 1;
//                                    if (newp1 == csdatabase.getP1_Position() && newp1prime == csdatabase.getP1prime_Position()) {
//                                        continue;
//                                    } else {
//                                        csdatabase.setP1_Position(newp1);
//                                        csdatabase.setP1prime_Position(newp1prime);
//                                        System.out.println(newp1);
//                                        System.out.println(newp1prime);
//                                        csdatabase.setCuration_Status("Cleavage site curated based on Uniprot protein sequence");
//                                        System.out.println("Cleavage site curated based on Uniprot protein sequence");
//                                    }
//
//                                } else {
//                                    csdatabase.setCuration_Status("Unmatched cleavage site; Cleavage site discarded");
//                                    System.out.println("Unmatched cleavage site; Cleavage site discarded");
//                                }
//
//                            }
//                            alldatanotcurated.add(csdatabase);
//                        }
//                    }
//                }
//            }
//        }
        //1b. Uniprot
//        String UniprotIDListURL = "http://www.uniprot.org/uniprot/?query=%22cleavage%3B%22+organism%3A9606+reviewed%3Ayes&format=list";
//        int lineCount = getLineCount(UniprotIDListURL);
//        String UniprotURLdebut = "http://www.uniprot.org/uniprot/?query=%22cleavage%3B%22+organism%3A9606+reviewed%3Ayes&format=xml&limit=100&offset=";
//
//        String UniprotURL = null;
//        int offset = 0;
//        while (offset < lineCount) {
//            UniprotURL = UniprotURLdebut + offset;
//            System.out.println(UniprotURL);
        String UniprotURL = "http://www.uniprot.org/uniprot/P05067.xml";

        NodeList entries = getEntries("/uniprot/entry[./feature[@type='site'][contains(@description, 'Cleavage')]]", parseUniprot(UniprotURL));
        for (int i = 0; i < entries.getLength(); i++) {

            SubstrateDatabaseEntry substratedatabase = new SubstrateDatabaseEntry();
            String uniprotid = getAccession(entries, i, substratedatabase);
            getUniSubstratepproteinname(entries, i, substratedatabase);
            String genename = getUniSubstrategenename(entries, i, substratedatabase);
            String substrateTaxon = getUniSubstratetaxonomy(entries, i, substratedatabase);
            String sequence = getUniSubstratesequence(entries, i, substratedatabase);
            String commentS = "-";

            NodeList csdbentries = getCsdbentries("./feature[@type='site'][contains(@description, 'Cleavage')]", entries.item(i));
            for (int j = 0; j < csdbentries.getLength(); j++) {

                Node n = csdbentries.item(j);
                LinkedList<String> descriptionlist = getInformation("./@description", n);
                String description = descriptionlist.getFirst();
                description = description.replaceAll("Cleavage; by ", "");
                description = description.replaceAll("Cleavage, first; by ", "");
                description = description.replaceAll("Cleavage, second; by ", "");
                description = description.replaceAll(",", ";");
                description = description.replaceAll("autolysis", genename);
                description = description.replaceAll("autocatalysis", genename);

                CsDatabaseEntry csdatabase = new CsDatabaseEntry();
                csdatabase.setSubstrate(substratedatabase);

                ProteaseDatabaseEntry proteasedatabase = new ProteaseDatabaseEntry();
//                proteasedatabase.setP_NL_Name(description);
                System.out.println(description);
                String proteaseTaxon = substrateTaxon;
                proteasedatabase.setP_Taxon(proteaseTaxon);
                String curationUni = "-";
                csdatabase.setCuration_Status(curationUni);

                String commentP = mapProteasetoLibrairy(commentS, proteaseTaxon, description, csdatabase, proteasedatabase);

                LinkedList<String> p1intlist = getInformation("./location/begin/@position", n);
                String p1int = null;
                int intP1 = 0;
                if (p1intlist.isEmpty()) {
                    LinkedList<String> positionlist = getInformation("./location//position/@position", n);
                    p1int = positionlist.getFirst();
                    System.out.println(p1int);
                    intP1 = Integer.parseInt(p1int);
                    csdatabase.setP1_Position(intP1);
                    char aaP1 = sequence.charAt(intP1 - 1);
                    String saaP1 = Character.toString(aaP1);
                    System.out.println(saaP1);
                    csdatabase.setP1_Sequence(saaP1);
                } else {
                    p1int = p1intlist.getFirst();
                    System.out.println(p1int);
                    intP1 = Integer.parseInt(p1int);
                    csdatabase.setP1_Position(intP1);
                    char aaP1 = sequence.charAt(intP1 - 1);
                    String saaP1 = Character.toString(aaP1);
                    System.out.println(aaP1);
                    csdatabase.setP1_Sequence(saaP1);
                }

                LinkedList<String> p1primeintlist = getInformation("./location/end/@position", n);
                String p1primeint = null;
                int intP1prime = 0;
                if (p1primeintlist.isEmpty()) {
                    intP1prime = intP1 + 1;
                    System.out.println(intP1prime);
                    csdatabase.setP1prime_Position(intP1prime);
                    char aaP1prime = sequence.charAt(intP1prime - 1);
                    String saaP1prime = Character.toString(aaP1prime);
                    System.out.println(aaP1prime);
                    csdatabase.setP1prime_Sequence(saaP1prime);

                } else {
                    p1primeint = p1primeintlist.getFirst();
                    System.out.println(p1primeint);
                    intP1prime = Integer.parseInt(p1primeint);
                    csdatabase.setP1prime_Position(intP1prime);
                    char aaP1prime = sequence.charAt(intP1prime - 1);
                    String saaP1prime = Character.toString(aaP1prime);
                    System.out.println(aaP1prime);
                    csdatabase.setP1prime_Sequence(saaP1prime);
                }

                int sequencelength = sequence.length();

                if (intP1 - 4 > 0 || intP1 - 4 == 0) {
                    if (intP1prime + 3 == sequencelength || intP1prime + 3 < sequencelength) {
                        String cleavagesite = sequence.substring(intP1 - 4, intP1prime + 3);
                        csdatabase.setCleavagesiteseaquence(cleavagesite);
                    } else {
                        String cleavagesite = sequence.substring(intP1 - 4, sequencelength);
                        csdatabase.setCleavagesiteseaquence(cleavagesite);
                    }
                } else {
                    if (intP1prime + 3 == sequencelength || intP1prime + 3 < sequencelength) {
                        String cleavagesite = sequence.substring(0, intP1prime + 3);
                        csdatabase.setCleavagesiteseaquence(cleavagesite);
                    } else {
                        String cleavagesite = sequence.substring(0, sequencelength);
                        csdatabase.setCleavagesiteseaquence(cleavagesite);
                    }
                }

                String csdburl = "http://www.uniprot.org/uniprot/" + uniprotid;
                csdatabase.setExternal_Link(csdburl);
                System.out.println(csdburl);
                SimpleDateFormat format = new SimpleDateFormat("dd-MM-yyyy");
                Calendar originalDate = Calendar.getInstance();
                String dateString = format.format(originalDate.getTime());
                System.out.println(dateString);
                csdatabase.setCreation_Date(dateString);
                csdatabase.setPMID("-");
                csdatabase.setComment(description);


                alldatanotcurated.add(csdatabase);
            }

//            }
//            offset = offset + 100;
        }

        getManualentries(alldatanotcurated);


        try {
            System.out.println("-----------------");
            csvWriter = new PrintStream("notcuratedProteasixDB" + "_" + version + ".csv");
//            populateHeaders(csvWriter);
            for (CsDatabaseEntry csdatabase : alldatanotcurated) {
                populateData(csvWriter, csdatabase);
            }

        } catch (FileNotFoundException ex) {
            Logger.getLogger(ProteasixCSdatabase.class.getName()).log(Level.SEVERE, null, ex);
        } finally {
            csvWriter.close();
        }


        BufferedReader bReader = createBufferedreader("//Users/julieklein/NetBeansProjects/ProteasixCSdatabase/notcuratedProteasixDB_AVRIL2012.csv");
        String line;

        Map<String, List<Set<String>>> hmap = new HashMap<String, List<Set<String>>>();
        while ((line = bReader.readLine()) != null) {
            String splitarray[] = line.split(",");
            String p_nl = splitarray[0];
            String p_symbol = splitarray[1];
            String p_uni = splitarray[2];
            String p_ec = splitarray[3];
            String p_taxon = splitarray[4];
            String s_nl = splitarray[5];
            String s_symbol = splitarray[6];
            String s_uni = splitarray[7];
            String s_taxon = splitarray[8];
            String cs = splitarray[9];
            String p1pos = splitarray[10];
            int ip1pos = Integer.parseInt(p1pos);
            String p1primepos = splitarray[11];
            int ip1primepos = Integer.parseInt(p1primepos);
            String p1seq = splitarray[12];
            String p1primeseq = splitarray[13];
            String ext_link = splitarray[14];
            String pmid = splitarray[15];
            String comment = splitarray[16];
            String curation = splitarray[17];
            String date = splitarray[18];
            String key = p_uni + s_uni + cs + p1pos;
            if (!hmap.containsKey(key) && !curation.contains("discarded")) {
                List value = new ArrayList<Set<String>>();
                for (int i = 0; i < 19; i++) {
                    value.add(new HashSet<String>());
                }
                hmap.put(key, value);
            }
            hmap.get(key).get(0).add(p_nl);
            hmap.get(key).get(1).add(p_symbol);
            hmap.get(key).get(2).add(p_uni);
            hmap.get(key).get(3).add(p_ec);
            hmap.get(key).get(4).add(p_taxon);
            hmap.get(key).get(5).add(s_nl);
            hmap.get(key).get(6).add(s_symbol);
            hmap.get(key).get(7).add(s_uni);
            hmap.get(key).get(8).add(s_taxon);
            hmap.get(key).get(9).add(cs);
            hmap.get(key).get(10).add(p1pos);
            hmap.get(key).get(11).add(p1primepos);
            hmap.get(key).get(12).add(p1seq);
            hmap.get(key).get(13).add(p1primeseq);
            hmap.get(key).get(14).add(ext_link);
            hmap.get(key).get(15).add(pmid);
            hmap.get(key).get(16).add(comment);
            hmap.get(key).get(17).add(curation);
            hmap.get(key).get(18).add(date);
        }


        Iterator iterator = hmap.values().iterator();
        while (iterator.hasNext()) {
            ProteaseDatabaseEntry proteasedatabase = new ProteaseDatabaseEntry();
            SubstrateDatabaseEntry substratedatabse = new SubstrateDatabaseEntry();
            CsDatabaseEntry csdatabasecur = new CsDatabaseEntry();
            csdatabasecur.setProtease(proteasedatabase);
            csdatabasecur.setSubstrate(substratedatabse);
            String values = iterator.next().toString();
//            System.out.println("-----------------");
//            System.out.println(values + "\n");

            String splitarray[] = values.split("\\], \\[");
            String p_nl = splitarray[0];
            p_nl = p_nl.replaceAll("\\[", "");
            String p_symbol = splitarray[1];
            String p_uni = splitarray[2];
            String p_ec = splitarray[3];
            String p_taxon = splitarray[4];
            String s_nl = splitarray[5];
            String s_symbol = splitarray[6];
            String s_uni = splitarray[7];
            String s_taxon = splitarray[8];
            String cs = splitarray[9];
            String p1pos = splitarray[10];
            int ip1pos = Integer.parseInt(p1pos);
            String p1primepos = splitarray[11];
            int ip1primepos = Integer.parseInt(p1primepos);
            String p1seq = splitarray[12];
            String p1primeseq = splitarray[13];
            String ext_link = splitarray[14];
            ext_link = ext_link.replaceAll(",", ";");
            String pmid = splitarray[15];
            pmid = pmid.replaceAll(",", ";");
            String comment = splitarray[16];
            comment = comment.replaceAll(",", ";");
            String curation = splitarray[17];
            curation.replaceAll(",", ";");
            String date = splitarray[18];
            date = date.replaceAll("\\]", "");
            proteasedatabase.setP_NL_Name(p_nl);
            proteasedatabase.setP_Symbol(p_symbol);
            proteasedatabase.setP_UniprotID(p_uni);
            proteasedatabase.setP_EC_Number(p_ec);
            proteasedatabase.setP_Taxon(p_taxon);
            substratedatabse.setS_NL_Name(s_nl);
            substratedatabse.setS_Symbol(s_symbol);
            substratedatabse.setS_UniprotID(s_uni);
            substratedatabse.setS_Taxon(s_taxon);
            csdatabasecur.setCleavagesiteseaquence(cs);
            csdatabasecur.setP1_Position(ip1pos);
            csdatabasecur.setP1prime_Position(ip1primepos);
            csdatabasecur.setP1_Sequence(p1seq);
            csdatabasecur.setP1prime_Sequence(p1primeseq);
            csdatabasecur.setCreation_Date(date);
            csdatabasecur.setCuration_Status(curation);
            csdatabasecur.setExternal_Link(ext_link);
            csdatabasecur.setPMID(pmid);
            csdatabasecur.setComment(comment);
            alldatacurated.add(csdatabasecur);
        }



//    
        try {
            System.out.println("-----------------");
            csvWriter = new PrintStream("curatedProteasixDB" + "_" + version + ".csv");
//            populateHeaders(csvWriter);
            for (CsDatabaseEntry csdatabasecur : alldatacurated) {
                populateData(csvWriter, csdatabasecur);
            }

        } catch (FileNotFoundException ex) {
            Logger.getLogger(ProteasixCSdatabase.class.getName()).log(Level.SEVERE, null, ex);
        } finally {
            csvWriter.close();
        }
//    }
//3. Populate in mysql with version
          try {
        String host = "jdbc:mysql://localhost:3306/";
        String dbName = "Proteasix" + version;
        String usermame = "root";
        String pwd = "kschoicesql";
        Connection conn = DriverManager.getConnection("jdbc:mysql://localhost/mysql?user=root&password=kschoicesql");
        Statement stat = conn.createStatement();
        stat.execute("CREATE DATABASE " + dbName);
        conn.close();
        
        Connection conn2 = DriverManager.getConnection(host + dbName + "?user=root&password=kschoicesql");
        Statement stat2 = conn2.createStatement();
        
        stat2.execute("CREATE TABLE PROTEASE(P_NL_Name VARCHAR(100), P_Symbol VARCHAR(100), P_UniprotID VARCHAR(100), P_EC_Number VARCHAR(100), P_Taxon VARCHAR(100))");
        stat2.execute("CREATE TABLE TAXONOMY(Taxon VARCHAR(100), PRIMARY KEY(Taxon))");
        stat2.execute("CREATE TABLE SUBSTRATE(S_NL_Name VARCHAR(100), S_Symbol VARCHAR(100), S_UniprotID VARCHAR(100), S_Taxon VARCHAR(100), PRIMARY KEY(S_NL_Name, S_Symbol, S_UniprotID, S_Taxon), FOREIGN KEY(S_Taxon) references TAXONOMY)");
        stat2.execute("CREATE TABLE CSDATABASE(P_NL_Name VARCHAR(100), P_Symbol VARCHAR(100), P_UniprotID VARCHAR(100), P_EC_Number VARCHAR(100), P_Taxon VARCHAR(100), S_NL_Name VARCHAR(100), S_Symbol VARCHAR(100), S_UniprotID VARCHAR(100), S_Taxon VARCHAR(100), Cleavage_Site VARCHAR(8), P1_Position INT, P1prime_Position INT, P1_Sequence VARCHAR(1), P1prime_Sequence VARCHAR(1), External_Link VARCHAR(1000), PMID VARCHAR(100), Comment VARCHAR(100), Curation_Status VARCHAR(100), Creation_Date VARCHAR(100))");
//        

        conn2.close();

        
        System.out.println("Connection Success");
          
        } catch (SQLException ex) {
            System.err.println("Connection Failed");
        }
    }

    private void populateHeaders(PrintStream csvWriter) throws FileNotFoundException {
        
        csvWriter.print("P_NL_Name");
        csvWriter.print(",");
        csvWriter.print("P_Symbol");
        csvWriter.print(",");
        csvWriter.print("P_UniprotID");
        csvWriter.print(",");
        csvWriter.print("P_EC_Number");
        csvWriter.print(",");
        csvWriter.print("P_Taxon");
        csvWriter.print(",");
        csvWriter.print("S_NL_Name");
        csvWriter.print(",");
        csvWriter.print("S_Symbol");
        csvWriter.print(",");
        csvWriter.print("S_UniprotID");
        csvWriter.print(",");
        csvWriter.print("S_Taxon");
//            csvWriter.print(",");
//            csvWriter.print("Substrate Sequence");
        csvWriter.print(",");
        csvWriter.print("Cleavage_Site");
        csvWriter.print(",");
        csvWriter.print("P1_Position");
        csvWriter.print(",");
        csvWriter.print("P1prime_Position");
        csvWriter.print(",");
        csvWriter.print("P1_Sequence");
        csvWriter.print(",");
        csvWriter.print("P1prime_Sequence");
        csvWriter.print(",");
        csvWriter.print("External_Link");
        csvWriter.print(",");
        csvWriter.print("PMID");
        csvWriter.print(",");
        csvWriter.print("Comment");
        csvWriter.print(",");
        csvWriter.print("Curation_Status");
        csvWriter.print(",");
        csvWriter.print("Creation_Date");
        csvWriter.print("\n");
    }

    private void populateData(PrintStream csvWriter, CsDatabaseEntry csdatabase) {
        //System.out.println(cleavageSiteDBEntry);
        
        csvWriter.print(csdatabase.protease.getP_NL_Name());
        csvWriter.print(",");
        csvWriter.print(csdatabase.protease.getP_Symbol());
        csvWriter.print(",");
        csvWriter.print(csdatabase.protease.getP_UniprotID());
        csvWriter.print(",");
        csvWriter.print(csdatabase.protease.getP_EC_Number());
        csvWriter.print(",");
        csvWriter.print(csdatabase.protease.getP_Taxon());
        csvWriter.print(",");
        csvWriter.print(csdatabase.substrate.getS_NL_Name());
        csvWriter.print(",");
        csvWriter.print(csdatabase.substrate.getS_Symbol());
        csvWriter.print(",");
        csvWriter.print(csdatabase.substrate.getS_UniprotID());
//                csvWriter.print(",");
//                csvWriter.print(csdatabase.substrate.getSubstratesequence());
        csvWriter.print(",");
        csvWriter.print(csdatabase.substrate.getS_Taxon());
        csvWriter.print(",");
        csvWriter.print(csdatabase.getCleavagesiteseaquence());
        csvWriter.print(",");
        csvWriter.print(csdatabase.getP1_Position());
        csvWriter.print(",");
        csvWriter.print(csdatabase.getP1prime_Position());
        csvWriter.print(",");
        csvWriter.print(csdatabase.getP1_Sequence());
        csvWriter.print(",");
        csvWriter.print(csdatabase.getP1prime_Sequence());
        csvWriter.print(",");
        csvWriter.print(csdatabase.getExternal_Link());
        csvWriter.print(",");
        csvWriter.print(csdatabase.getPMID());
        csvWriter.print(",");
        csvWriter.print(csdatabase.getComment());
        csvWriter.print(",");
        csvWriter.print(csdatabase.getCuration_Status());
        csvWriter.print(",");
        csvWriter.print(csdatabase.getCreation_Date());
        csvWriter.print("\n");
    }

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws MalformedURLException, IOException {
        // TODO code application logic here
        ProteasixCSdatabase ProteasixCSdatabase = new ProteasixCSdatabase();
    }
}
