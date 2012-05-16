/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package proteasixcsdatabase;

/**
 *
 * @author julieklein
 */
public class SubstrateDatabaseEntry {

    public String getS_NL_Name() {
        return S_NL_Name;
    }

    public void setS_NL_Name(String S_NL_Name) {
        this.S_NL_Name = S_NL_Name;
    }

    public String getS_Symbol() {
        return S_Symbol;
    }

    public void setS_Symbol(String S_Symbol) {
        this.S_Symbol = S_Symbol;
    }

    public String getS_Taxon() {
        return S_Taxon;
    }

    public void setS_Taxon(String S_Taxon) {
        this.S_Taxon = S_Taxon;
    }

    public String getS_UniprotID() {
        return S_UniprotID;
    }

    public void setS_UniprotID(String S_UniprotID) {
        this.S_UniprotID = S_UniprotID;
    }
       
    private String S_NL_Name;
    private String S_Symbol;
    private String S_UniprotID;
    private String S_Taxon;
    
    
    
}
