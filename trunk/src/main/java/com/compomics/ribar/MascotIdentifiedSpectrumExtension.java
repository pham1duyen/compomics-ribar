package com.compomics.ribar;

import com.compomics.mslims.util.mascot.MascotIdentifiedSpectrum;

import java.util.Vector;

/**
 * Created by IntelliJ IDEA.
 * User: Niklaas
 * Date: 29/04/11
 * Time: 14:13
 * To change this template use File | Settings | File Templates.
 */
public class MascotIdentifiedSpectrumExtension extends MascotIdentifiedSpectrum {

    double iMsMsSpectrumIntensity;


    public MascotIdentifiedSpectrumExtension(double aMsMsSpectrumIntensity) {
        iMsMsSpectrumIntensity = aMsMsSpectrumIntensity;
    }

    public double getTotalIntensity() {
        return iMsMsSpectrumIntensity;
    }

    public double getMatchedIntensity() {
        //ToDo not implemented yet
        return 0.0;
    }

    public Vector<String> getAccessions(){
        Vector<String> lLocalAccessions = new Vector<String>();
        String lAccession = this.getAccession(null);
        lLocalAccessions.add(lAccession);
        // See if there are any isoforms in the description.
        String tempDesc = this.getDescription(lAccession);
        String isoforms = this.getIsoformAccessions(lAccession);

        int startIsoforms = -1;
        boolean lMoreFound = false;
        if ((startIsoforms = tempDesc.indexOf("^A")) >= 0) {
            lMoreFound = true;
            String tempDesc2 = tempDesc.substring(0, startIsoforms);
            lLocalAccessions.add(tempDesc2.substring(0, tempDesc2.indexOf(" ")));
            this.setDescription(tempDesc2, lAccession);
            if (isoforms == null) {
                isoforms = tempDesc.substring(startIsoforms + 2);
            } else {
                isoforms += tempDesc.substring(startIsoforms);
            }
        }
        if(!lMoreFound){
            if(isoforms.length()>0){
                lLocalAccessions.add(isoforms.substring(0,isoforms.indexOf(" ")));
            }
        }
        return lLocalAccessions;

    }

}
