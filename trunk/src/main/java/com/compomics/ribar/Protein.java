package com.compomics.ribar;

import com.compomics.mslims.util.mascot.MascotIdentifiedSpectrum;
import com.compomics.rover.general.sequenceretriever.UniprotSequenceRetriever;
import com.compomics.util.protein.AASequenceImpl;
import com.compomics.util.protein.Enzyme;
import com.compomics.util.protein.Header;
import org.apache.commons.math.stat.descriptive.DescriptiveStatistics;

import java.util.Vector;

/**
 * Created by IntelliJ IDEA.
 * User: Niklaas
 * Date: 29/04/11
 * Time: 14:33
 * To change this template use File | Settings | File Templates.
 */
public class Protein {

    private String iAccession;
    //private UniprotProtein iUniprotProtein;
    private Vector<MascotIdentifiedSpectrumExtension> iIdentifications = new Vector<MascotIdentifiedSpectrumExtension>();
    private com.compomics.util.protein.Protein[] iPeptides;
    private boolean[] iSelectedPeptides;
    private int iLength = 0;
    private int iNumberOfSequences = 0;
    private int iNumberOfUniqueParentIons = 0;
    private int iNumberOfIdentifiablePeptides = 0;
    private Project iParent;
    private com.compomics.util.protein.Protein iUtilProtein;
    private double iMass = 0.0;
    private String iSequence = null;
    private Vector<Vector<MascotIdentifiedSpectrumExtension>> iGroupedIdentificationsBySequence = new Vector<Vector<MascotIdentifiedSpectrumExtension>>();
    private Vector<Vector<MascotIdentifiedSpectrumExtension>> iGroupedIdentificationsByModifiedSequence = new Vector<Vector<MascotIdentifiedSpectrumExtension>>();
    private double iNSAF = 0.0;


    public Protein(String lAccession, Project aProject){
        this.iAccession = lAccession;
        this.iParent = aProject;
    }


    public void matchIdentification(Vector<MascotIdentifiedSpectrumExtension> lIdentifications){
        for(int i = 0; i<lIdentifications.size(); i ++){
            Vector<String> lAccessions = lIdentifications.get(i).getAccessions();
            for(int j = 0; j<lAccessions.size(); j ++){
                if(lAccessions.get(j).indexOf(iAccession) >= 0){
                    iIdentifications.add(lIdentifications.get(i));
                }
            }
            /*else if(lIdentifications.get(i).getIsoforms().indexOf(iAccession) >= 0){
                if(iAccession.indexOf("rev_") == -1 ){
                    if(lIdentifications.get(i).getIsoforms().indexOf("rev_" + iAccession) == -1){
                        iIdentifications.add(lIdentifications.get(i));
                    }
                } else {
                    iIdentifications.add(lIdentifications.get(i));
                }
            }   */
        }

        countNumberOfSequences();
        countNumberOfUniqueParentIons();
    }

    public void selectIdentifiablePeptides(double lLower, double lUpper){
        if(iSequence == null){
            String lSequence = null;
            try {
                UniprotSequenceRetriever lUn = new UniprotSequenceRetriever(iAccession);
                lSequence = lUn.getSequence();

            } catch (Exception e) {
                e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            }
            this.iUtilProtein = new com.compomics.util.protein.Protein(Header.parseFromFASTA(""), new AASequenceImpl(lSequence));
            iMass = iUtilProtein.getSequence().getMass();
            Enzyme lEnzyme = new Enzyme("Trypsin", "KR", "P", "Cterm", 0);
            iPeptides =  lEnzyme.cleave(iUtilProtein);
            iSelectedPeptides =  new boolean[iPeptides.length];
            iLength = lSequence.length();
            iSequence = lSequence;
        }

        int lCounter = 0;
        for(int i = 0; i<iPeptides.length; i ++){
            if(lLower <= iPeptides[i].getSequence().getMass() && iPeptides[i].getSequence().getMass() <= lUpper){
                iSelectedPeptides[i] = true;
                lCounter = lCounter + 1;
            } else {
                iSelectedPeptides[i] = false;
            }
        }
        iNumberOfIdentifiablePeptides = lCounter;
    }


    public int countNumberOfSequences(){
        if(iNumberOfSequences == 0){
            Vector<String> lSequences = new Vector<String>();
            for(int i = 0; i<iIdentifications.size(); i ++){
                boolean lFound = false;

                for(int j = 0; j<lSequences.size(); j ++){
                    if(iIdentifications.get(i).getSequence().equalsIgnoreCase(lSequences.get(j))){
                        lFound = true;
                        iGroupedIdentificationsBySequence.get(j).add(iIdentifications.get(i));
                        j = lSequences.size();
                    }
                }

                if(!lFound){
                    lSequences.add(iIdentifications.get(i).getSequence());
                    Vector<MascotIdentifiedSpectrumExtension> lTemp = new Vector<MascotIdentifiedSpectrumExtension>();
                    lTemp.add(iIdentifications.get(i));
                    iGroupedIdentificationsBySequence.add(lTemp);
                }
            }
            iNumberOfSequences = lSequences.size();
        }
        return iNumberOfSequences;
    }


    public int countNumberOfUniqueParentIons(){
        if(iNumberOfUniqueParentIons == 0){
            Vector<String> lSequences = new Vector<String>();
            for(int i = 0; i<iIdentifications.size(); i ++){
                boolean lFound = false;

                for(int j = 0; j<lSequences.size(); j ++){
                    //String lTempSeq = iIdentifications.get(i).getModified_sequence() + "_" + iIdentifications.get(i).getCharge();
                    String lTempSeq = iIdentifications.get(i).getModifiedSequence();
                    if(lTempSeq.equalsIgnoreCase(lSequences.get(j))){
                        lFound = true;
                        iGroupedIdentificationsByModifiedSequence.get(j).add(iIdentifications.get(i));
                        j = lSequences.size();
                    }
                }

                if(!lFound){
                    //lSequences.add(iIdentifications.get(i).getModified_sequence() + "_" + iIdentifications.get(i).getCharge());
                    lSequences.add(iIdentifications.get(i).getModifiedSequence());
                    Vector<MascotIdentifiedSpectrumExtension> lTemp = new Vector<MascotIdentifiedSpectrumExtension>();
                    lTemp.add(iIdentifications.get(i));
                    iGroupedIdentificationsByModifiedSequence.add(lTemp);
                }
            }
            iNumberOfUniqueParentIons = lSequences.size();
        }
        return iNumberOfUniqueParentIons;
    }


    public double getPAI(boolean lUniqueParentIons){
        double lPAI;
        if(lUniqueParentIons){
            lPAI = (double) iNumberOfUniqueParentIons / (double) iNumberOfIdentifiablePeptides;
        } else {
            lPAI = (double) iNumberOfSequences / (double) iNumberOfIdentifiablePeptides;
        }
        if(lPAI >= 1.0){
            //System.out.println("Trouble?");
        }
        return lPAI;
    }

    public double getEmPAI(boolean lUniqueParentIons){
        double lPAI = this.getPAI(lUniqueParentIons);
        double lEmPAI = Math.pow(10, lPAI) - 1.0;
        return lEmPAI;
    }


    public double getProteinContentInPercMol(boolean lUniqueParentIons){
        double lEmPAISum = iParent.getEmPAISum(lUniqueParentIons);
        double lResult = (this.getEmPAI(lUniqueParentIons) / lEmPAISum) * 100.0;
        return lResult;
    }

    public double getMass(){
        return iMass;
    }

    public double getProteinContentInPercWeight(boolean lModifiedSequences){
        double lEmPAIWeightSum = iParent.getEmPAIWeightSum(lModifiedSequences);
        double lResult = ( (this.getEmPAI(lModifiedSequences) * this.getMass()) / lEmPAIWeightSum) * 100.0;
        return lResult;
    }

    public double getSI(boolean lTotal){

        double lSum = 0.0;
        for(int i = 0; i<iIdentifications.size(); i ++){
            if(lTotal){
                lSum = lSum + iIdentifications.get(i).getTotalIntensity();
            } else {
                lSum = lSum + iIdentifications.get(i).getMatchedIntensity();
            }
        }
        return lSum;
    }

    public double getSIForModPeptide(boolean lTotal, String lModPeptide){

        double lSum = 0.0;
        for(int i = 0; i<iIdentifications.size(); i ++){
            if(lModPeptide.equalsIgnoreCase(iIdentifications.get(i).getModifiedSequence())){
                if(lTotal){
                    lSum = lSum + iIdentifications.get(i).getTotalIntensity();
                } else {
                    lSum = lSum + iIdentifications.get(i).getMatchedIntensity();
                }
            }
        }
        return lSum;
    }

    public double getSIN(boolean lTotal){

        double lSiSum = iParent.getSiSum(lTotal);
        double lSin = (this.getSI(lTotal)/lSiSum) / (double) iLength;
        return lSin;
    }



    public double getSINPAI(boolean lTotal){

        double lSiSum = iParent.getSiSum(lTotal);
        double lSin = ((this.getSI(lTotal)*this.getPAI(true))/lSiSum) / (double) iLength;
        return lSin;

    }

    public double getSPL(){
        return (double) iIdentifications.size() / (double) iLength;
    }

    public int countNumberOfSpectra(){
        return iIdentifications.size();
    }

    public double getNSAF(){
        if(iNSAF == 0.0){
            double lSPLSum = iParent.getSPLSum();
            double lNSAF = this.getSPL() / lSPLSum;
            iNSAF = lNSAF;
        }
        return iNSAF;
    }

    public double getPeptideSpectrumMatches(){
        return (double) iIdentifications.size();
    }


    public String getAccession() {
        return iAccession;
    }

    public int getLength(){
        return iLength;
    }

    public double getSINQualityMean(boolean lModifiedSequences){
        DescriptiveStatistics lNorm = new DescriptiveStatistics();
        if(lModifiedSequences){
            for(int i = 0; i<iGroupedIdentificationsByModifiedSequence.size(); i ++){
                Vector<MascotIdentifiedSpectrumExtension> lIds = iGroupedIdentificationsByModifiedSequence.get(i);
                DescriptiveStatistics lTemp = new DescriptiveStatistics();
                for(int j = 0; j<lIds.size(); j ++){
                    lTemp.addValue(lIds.get(j).getTotalIntensity());
                }
                if(lIds.size() > 0){
                    lNorm.addValue(lTemp.getStandardDeviation() / lTemp.getMean());
                } else {
                    lNorm.addValue(0.0);
                }
            }
        } else {
            for(int i = 0; i<iGroupedIdentificationsBySequence.size(); i ++){
                Vector<MascotIdentifiedSpectrumExtension> lIds = iGroupedIdentificationsBySequence.get(i);
                DescriptiveStatistics lTemp = new DescriptiveStatistics();
                for(int j = 0; j<lIds.size(); j ++){
                    lTemp.addValue(lIds.get(j).getTotalIntensity());
                }
                if(lIds.size() > 0){
                    lNorm.addValue(lTemp.getStandardDeviation() / lTemp.getMean());
                } else {
                    lNorm.addValue(0.0);
                }
            }
        }
        return lNorm.getMean();
    }


    public DescriptiveStatistics getStatsForPeptide(boolean lModifiedSequences, String aPeptide){
        DescriptiveStatistics lNorm = new DescriptiveStatistics();
        if(lModifiedSequences){
            for(int i = 0; i<iGroupedIdentificationsByModifiedSequence.size(); i ++){
                if(iGroupedIdentificationsByModifiedSequence.get(i).get(0).getModifiedSequence().equalsIgnoreCase(aPeptide)){
                    Vector<MascotIdentifiedSpectrumExtension> lIds = iGroupedIdentificationsByModifiedSequence.get(i);
                    for(int j = 0; j<lIds.size(); j ++){
                        lNorm.addValue(lIds.get(j).getTotalIntensity());
                    }
                }
            }
        } else {
            for(int i = 0; i<iGroupedIdentificationsBySequence.size(); i ++){
                if(iGroupedIdentificationsBySequence.get(i).get(0).getModifiedSequence().equalsIgnoreCase(aPeptide)){
                    Vector<MascotIdentifiedSpectrumExtension> lIds = iGroupedIdentificationsBySequence.get(i);
                    for(int j = 0; j<lIds.size(); j ++){
                        lNorm.addValue(lIds.get(j).getTotalIntensity());
                    }
                }
            }
        }
        return lNorm;
    }

    public double getSINQualitySD(boolean lModifiedSequences){
        DescriptiveStatistics lNorm = new DescriptiveStatistics();
        if(lModifiedSequences){
            for(int i = 0; i<iGroupedIdentificationsByModifiedSequence.size(); i ++){
                Vector<MascotIdentifiedSpectrumExtension> lIds = iGroupedIdentificationsByModifiedSequence.get(i);
                DescriptiveStatistics lTemp = new DescriptiveStatistics();
                for(int j = 0; j<lIds.size(); j ++){
                    lTemp.addValue(lIds.get(j).getTotalIntensity());
                }
                if(lIds.size() > 0){
                    lNorm.addValue(lTemp.getStandardDeviation() / lTemp.getMean());
                } else {
                    lNorm.addValue(0.0);
                }
            }
        } else {
            for(int i = 0; i<iGroupedIdentificationsBySequence.size(); i ++){
                Vector<MascotIdentifiedSpectrumExtension> lIds = iGroupedIdentificationsBySequence.get(i);
                DescriptiveStatistics lTemp = new DescriptiveStatistics();
                for(int j = 0; j<lIds.size(); j ++){
                    lTemp.addValue(lIds.get(j).getTotalIntensity());
                }
                if(lIds.size() > 0){
                    lNorm.addValue(lTemp.getStandardDeviation() / lTemp.getMean());
                } else {
                    lNorm.addValue(0.0);
                }
            }
        }
        return lNorm.getStandardDeviation();
    }

    public double getCVofMeanIntensities(boolean lModifiedSequences){
        DescriptiveStatistics lNorm = new DescriptiveStatistics();
        if(lModifiedSequences){
            for(int i = 0; i<iGroupedIdentificationsByModifiedSequence.size(); i ++){
                Vector<MascotIdentifiedSpectrumExtension> lIds = iGroupedIdentificationsByModifiedSequence.get(i);
                DescriptiveStatistics lTemp = new DescriptiveStatistics();
                for(int j = 0; j<lIds.size(); j ++){
                    lTemp.addValue(lIds.get(j).getTotalIntensity());
                }
                if(lIds.size() > 0){
                    lNorm.addValue(lTemp.getMean());
                }
            }
        } else {
            for(int i = 0; i<iGroupedIdentificationsBySequence.size(); i ++){
                Vector<MascotIdentifiedSpectrumExtension> lIds = iGroupedIdentificationsBySequence.get(i);
                DescriptiveStatistics lTemp = new DescriptiveStatistics();
                for(int j = 0; j<lIds.size(); j ++){
                    lTemp.addValue(lIds.get(j).getTotalIntensity());
                }
                if(lIds.size() > 0){
                    lNorm.addValue(lTemp.getMean());
                }
            }
        }
        if(lNorm.getValues().length <= 1){
            return Double.NaN;
        }
        return lNorm.getStandardDeviation()/lNorm.getMean();
    }

    public double getSDofIntensities(boolean lModifiedSequences){
        DescriptiveStatistics lNorm = new DescriptiveStatistics();
        if(lModifiedSequences){
            for(int i = 0; i<iGroupedIdentificationsByModifiedSequence.size(); i ++){
                Vector<MascotIdentifiedSpectrumExtension> lIds = iGroupedIdentificationsByModifiedSequence.get(i);
                for(int j = 0; j<lIds.size(); j ++){
                    lNorm.addValue(lIds.get(j).getTotalIntensity());
                }
            }
        } else {
            for(int i = 0; i<iGroupedIdentificationsBySequence.size(); i ++){
                Vector<MascotIdentifiedSpectrumExtension> lIds = iGroupedIdentificationsBySequence.get(i);
                for(int j = 0; j<lIds.size(); j ++){
                    lNorm.addValue(lIds.get(j).getTotalIntensity());
                }
            }
        }
        if(lNorm.getValues().length <= 1){
            return Double.NaN;
        }
        return lNorm.getStandardDeviation();
    }

    public double getSDofMeanIntensities(boolean lModifiedSequences){
        DescriptiveStatistics lNorm = new DescriptiveStatistics();
        if(lModifiedSequences){
            for(int i = 0; i<iGroupedIdentificationsByModifiedSequence.size(); i ++){
                Vector<MascotIdentifiedSpectrumExtension> lIds = iGroupedIdentificationsByModifiedSequence.get(i);
                DescriptiveStatistics lTemp = new DescriptiveStatistics();
                for(int j = 0; j<lIds.size(); j ++){
                    lTemp.addValue(lIds.get(j).getTotalIntensity());
                }
                if(lIds.size() > 0){
                    lNorm.addValue(lTemp.getMean());
                }
            }
        } else {
            for(int i = 0; i<iGroupedIdentificationsBySequence.size(); i ++){
                Vector<MascotIdentifiedSpectrumExtension> lIds = iGroupedIdentificationsBySequence.get(i);
                DescriptiveStatistics lTemp = new DescriptiveStatistics();
                for(int j = 0; j<lIds.size(); j ++){
                    lTemp.addValue(lIds.get(j).getTotalIntensity());
                }
                if(lIds.size() > 0){
                    lNorm.addValue(lTemp.getMean());
                }
            }
        }
        if(lNorm.getValues().length <= 1){
            return Double.NaN;
        }
        return lNorm.getStandardDeviation();
    }


    public void setSequence(String lSequence) {
        this.iSequence = lSequence;
        this.iUtilProtein = new com.compomics.util.protein.Protein(Header.parseFromFASTA(""), new AASequenceImpl(lSequence));
        iMass = iUtilProtein.getSequence().getMass();
        Enzyme lEnzyme = new Enzyme("Trypsin", "KR", "P", "Cterm", 0);
        iPeptides =  lEnzyme.cleave(iUtilProtein);
        iSelectedPeptides =  new boolean[iPeptides.length];
        iLength = iSequence.length();
    }

    public Vector<Vector<MascotIdentifiedSpectrumExtension>> getGroupedIdentificationsBySequence(){
        return iGroupedIdentificationsBySequence;
    }

    public Vector<Vector<MascotIdentifiedSpectrumExtension>> getGroupedIdentificationsByModifiedSequence(){
        return iGroupedIdentificationsByModifiedSequence;
    }


    public int getNumberOfObservablePeptides(){
        return iNumberOfIdentifiablePeptides;
    }

    public double getSSSC() {
        return this.getPeptideSpectrumMatches() / iParent.getPeptideSpectrumMatchesSum();
    }
}

