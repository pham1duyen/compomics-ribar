package com.compomics.ribar;

import com.compomics.mslims.util.mascot.MascotIdentifiedSpectrum;

import java.io.FileReader;
import java.io.IOException;
import java.io.LineNumberReader;
import java.sql.SQLException;
import java.util.Arrays;
import java.util.Collections;
import java.util.Vector;

/**
 * Created by IntelliJ IDEA.
 * User: Niklaas
 * Date: 29/04/11
 * Time: 14:31
 * To change this template use File | Settings | File Templates.
 */
public class Project {

    private Vector<Protein> iProteins = new Vector<Protein>();
    private double iEmPAISum = 0.0;
    private double iEmPAIWeightSum = 0.0;
    private double iSiSum = 0.0;
    private double iSPLSum = 0.0;


    public Project(Vector<MascotIdentifiedSpectrumExtension> aIdentifiedSpectrumExtensions) {

        //Get the protein accessions from all the identifications
        Vector<String> lAccessions = new Vector<String>();
        Vector<String> tempAccessions = new Vector<String>();
        // First cycle all to retrieve accession numbers from unique identifications.
        int liSize = aIdentifiedSpectrumExtensions.size();
        for (int i = 0; i < liSize; i++) {
            // Get the identified spectrum.
            MascotIdentifiedSpectrum lSpectrum = aIdentifiedSpectrumExtensions.elementAt(i);
            if (lSpectrum.getIsoformCount() == 1) {
                String lAccession = lSpectrum.getAccession(null);
                if (!tempAccessions.contains(lAccession)) {
                    tempAccessions.add(lAccession);
                }
            }
        }

        // Store all retrieved accession numbers.
        String[] accessions = new String[tempAccessions.size()];
        tempAccessions.toArray(accessions);
        Arrays.sort(accessions);
        // Invert the array.
        String[] accessionsInv = new String[accessions.length];
        for (int i = 0; i < accessions.length; i++) {
            accessionsInv[i] = accessions[accessions.length - (i + 1)];
        }

        for (int i = 0; i < aIdentifiedSpectrumExtensions.size(); i++) {

            Vector<String> lLocalAccessions = aIdentifiedSpectrumExtensions.get(i).getAccessions();


            for (int j = 0; j < lLocalAccessions.size(); j++) {
                boolean lFound = false;
                for(int k = 0; k<lAccessions.size(); k++){

                    if (lLocalAccessions.get(j).equalsIgnoreCase(lAccessions.get(k))) {
                        lFound = true;
                        k = lAccessions.size();
                    }
                }

                if (!lFound) {
                    lAccessions.add(lLocalAccessions.get(j));
                }
            }

            /*if(iIdentifications.get(i).getIsoforms().length() > 0){
                String[] lIsoforms = iIdentifications.get(i).getIsoforms().split("^A");

                for(int k = 0; k<lIsoforms.length ; k ++){
                    if(lIsoforms[k].length() > 0){

                        lIsoforms[k] = lIsoforms[k].substring(0, lIsoforms[k].indexOf(" "));
                        if(lIsoforms[k].length() == 5){
                            lIsoforms[k] = "A" + lIsoforms[k];
                        }
                        lFound = false;

                        for(int j = 0; j<lAccessions.size(); j ++){
                            if(lIsoforms[k].equalsIgnoreCase(lAccessions.get(j))){
                                lFound = true;
                                j = lAccessions.size();
                            }
                        }

                        if(!lFound){
                            lAccessions.add(lIsoforms[k]);
                        }
                    }
                }
            }*/
        }

        //Create all the proteins
        for (int i = 0; i < lAccessions.size(); i++) {
            Protein lProtein = new Protein(lAccessions.get(i), this);
            lProtein.matchIdentification(aIdentifiedSpectrumExtensions);
            iProteins.add(lProtein);
        }

    }


    public double getEmPAISum(boolean lUniqueParentIons) {
        if (iEmPAISum == 0.0) {
            double lSum = 0.0;
            for (int i = 0; i < iProteins.size(); i++) {
                double lEmpai = iProteins.get(i).getEmPAI(lUniqueParentIons);
                if (!Double.isNaN(lEmpai)) {
                    lSum = lSum + lEmpai;
                }
            }
            iEmPAISum = lSum;
        }
        return iEmPAISum;
    }

    public double getEmPAIWeightSum(boolean lUniqueParentIons) {
        if (iEmPAIWeightSum == 0.0) {
            double lSum = 0.0;
            for (int i = 0; i < iProteins.size(); i++) {
                double lEmpai = iProteins.get(i).getEmPAI(lUniqueParentIons) * iProteins.get(i).getMass();
                if (!Double.isNaN(lEmpai)) {
                    lSum = lSum + lEmpai;
                }
            }
            iEmPAIWeightSum = lSum;
        }
        return iEmPAIWeightSum;
    }

    public double getSiSum(boolean lTotal) {

        double lSum = 0.0;
        for (int i = 0; i < iProteins.size(); i++) {
            lSum = lSum + (iProteins.get(i).getSI(lTotal));
        }
        return lSum;

    }

    public double getPeptideSpectrumMatchesSum() {

        double lSum = 0.0;
        for (int i = 0; i < iProteins.size(); i++) {
            lSum = lSum + (iProteins.get(i).getPeptideSpectrumMatches());
        }
        return lSum;

    }

    public Vector<Protein> getProteins() {
        return iProteins;
    }

    public double getSPLSum() {
        if (iSPLSum == 0.0) {
            double lSum = 0.0;
            for (int i = 0; i < iProteins.size(); i++) {
                double lSPL = (iProteins.get(i).getSPL());
                if (!Double.isNaN(lSPL)) {
                    lSum = lSum + lSPL;
                }
            }
            iSPLSum = lSum;
        }
        return iSPLSum;
    }
}