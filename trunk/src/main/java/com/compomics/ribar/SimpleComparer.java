package com.compomics.ribar;

import com.compomics.mslims.util.mascot.MascotIdentifiedSpectrum;
import org.apache.commons.math.stat.descriptive.DescriptiveStatistics;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.Vector;

/**
 * Created by IntelliJ IDEA.
 * User: Niklaas
 * Date: 29/04/11
 * Time: 16:40
 * To change this template use File | Settings | File Templates.
 */
public class SimpleComparer {



    public SimpleComparer(Project lProject1, Project lProject2, String location){

        Vector<Protein> lProteins1 = lProject1.getProteins();
        Vector<Protein> lProteins2 = lProject2.getProteins();


        //create the writer
        try{


            BufferedWriter out = new BufferedWriter(new FileWriter(location));
            out.write("Protein,overlap,NSpectra1,Nions1,NSpectra2,Nions2,EmPAI ratio,SInMatch,NSAF ratio,SC,RIBAR,CoRIBAR,upsStatus,rev,extraRibar,notused\n");


            int lCounter = 0;
            for(int i = 0; i<lProteins1.size(); i ++){

                for(int j = 0; j<lProteins2.size(); j ++){
                    if(lProteins1.get(i).getAccession().equalsIgnoreCase(lProteins2.get(j).getAccession())){

                        Protein lP1 = lProteins1.get(i);
                        Protein lP2 = lProteins2.get(j);
                        String lUPSstatus = "0";
                        String lRevString = "0";
                        if(lProteins1.get(i).getAccession().endsWith("ups")){
                            lUPSstatus = "1";
                        }
                        if(lProteins1.get(i).getAccession().indexOf("rev") >= 0){
                            lRevString = "1";
                        }
                        Vector<Vector<MascotIdentifiedSpectrumExtension>> lIds1 = lProteins1.get(i).getGroupedIdentificationsByModifiedSequence();
                        Vector<Vector<MascotIdentifiedSpectrumExtension>> lIds2 = lProteins2.get(j).getGroupedIdentificationsByModifiedSequence();

                        Vector<DescriptiveStatistics> lStats1 = new Vector<DescriptiveStatistics>();
                        Vector<DescriptiveStatistics> lStats2 = new Vector<DescriptiveStatistics>();
                        Vector<DescriptiveStatistics> lNotUsedStats1 = new Vector<DescriptiveStatistics>();
                        Vector<DescriptiveStatistics> lNotUsedStats2 = new Vector<DescriptiveStatistics>();
                        Vector<String> lPeptides = new Vector<String>();
                        for(int k = 0; k<lIds1.size(); k ++){
                            for(int l = 0; l<lIds2.size(); l ++){
                                if(lIds1.get(k).get(0).getModifiedSequence().equalsIgnoreCase(lIds2.get(l).get(0).getModifiedSequence())){
                                    lPeptides.add(lIds1.get(k).get(0).getModifiedSequence());
                                    lStats1.add(lProteins1.get(i).getStatsForPeptide(true, lIds1.get(k).get(0).getModifiedSequence()));
                                    lStats2.add(lProteins2.get(j).getStatsForPeptide(true, lIds2.get(l).get(0).getModifiedSequence()));
                                }
                            }
                        }

                        double lIntSum1NotUsedInRibar = 0.0;
                        double lIntSum2NotUsedInRibar = 0.0;
                        double lCount1NotUsedInRibar = 0.0;
                        double lCount2NotUsedInRibar = 0.0;

                        for(int k = 0; k<lIds1.size(); k ++){
                            boolean lFound = false;
                            for(int l = 0; l<lIds2.size(); l ++){
                                if(lIds1.get(k).get(0).getModifiedSequence().equalsIgnoreCase(lIds2.get(l).get(0).getModifiedSequence())){
                                    lFound = true;
                                }
                            }
                            if(!lFound){
                                lIntSum1NotUsedInRibar = lIntSum1NotUsedInRibar + lProteins1.get(i).getStatsForPeptide(true, lIds1.get(k).get(0).getModifiedSequence()).getSum();
                                lCount1NotUsedInRibar = lCount1NotUsedInRibar + lProteins1.get(i).getStatsForPeptide(true, lIds1.get(k).get(0).getModifiedSequence()).getValues().length;
                            }
                        }
                        for(int k = 0; k<lIds2.size(); k ++){
                            boolean lFound = false;
                            for(int l = 0; l<lIds1.size(); l ++){
                                if(lIds2.get(k).get(0).getModifiedSequence().equalsIgnoreCase(lIds1.get(l).get(0).getModifiedSequence())){
                                    lFound = true;
                                }
                            }
                            if(!lFound){
                                lIntSum2NotUsedInRibar = lIntSum2NotUsedInRibar + lProteins2.get(j).getStatsForPeptide(true, lIds2.get(k).get(0).getModifiedSequence()).getSum();
                                lCount2NotUsedInRibar = lCount2NotUsedInRibar + lProteins2.get(j).getStatsForPeptide(true, lIds2.get(k).get(0).getModifiedSequence()).getValues().length;
                            }
                        }


                        DescriptiveStatistics lRibarStat = new DescriptiveStatistics();
                        DescriptiveStatistics lCoRibarStat = new DescriptiveStatistics();
                        double lIntens1 = 0.0;
                        double lIntens2 = 0.0;
                        int lCounter1 = 0;
                        int lCounter2 = 0;

                        for(int k = 0; k<lStats1.size(); k ++){

                            lIntens1 = lIntens1 + lStats1.get(k).getSum();
                            lIntens2 = lIntens2 + lStats2.get(k).getSum();
                            lCounter1 = lCounter1 + (int) lStats1.get(k).getN();
                            lCounter2 = lCounter2 + (int) lStats2.get(k).getN();

                            lRibarStat.addValue(Math.log(lStats1.get(k).getSum() / lStats2.get(k).getSum())/ Math.log(2));
                            lCoRibarStat.addValue(Math.log((lStats1.get(k).getSum()) / (lStats2.get(k).getSum()))/ Math.log(2));
                            //lRibarStat.addValue(Math.log(lStats1.get(k).getMean() / lStats2.get(k).getMean())/ Math.log(2));
                            //lRatioStatLog2Median.addValue(Math.log(BasicStats.median(lStats1.get(k).getValues(), false) / BasicStats.median(lStats2.get(k).getValues(), false)) / Math.log(2));
                            //lRibarStat.addValue( Math.log(   (BasicStats.median(lStats1.get(k).getValues(), false) * (double) lStats1.get(k).getN())  / (BasicStats.median(lStats2.get(k).getValues(), false)* (double) lStats2.get(k).getN()) ) /  Math.log(2));

                        }

                        //add another ratio with the rest of the intensities

                        /*
                        if(lCounter1 >= 1 || lCounter2 >= 1){

                            double lAverageInt1;
                            double lAverageInt2;
                            if(lCounter1 >= 1){
                                lAverageInt1 = lIntens1 / (double) lCounter1;
                                int lNumUnusedSpectra1 = lProteins1.get(i).countNumberOfSpectra() - lCounter1 + 1;
                                lAverageInt1 = lAverageInt1 * (double) lNumUnusedSpectra1;
                            } else {
                                lAverageInt1 = 1.0;
                            }

                            if(lCounter2 >= 1){
                                lAverageInt2 = lIntens2 / (double) lCounter2;
                                int lNumUnusedSpectra2 = lProteins2.get(j).countNumberOfSpectra() - lCounter2 + 1;
                                lAverageInt2 = lAverageInt2 * (double) lNumUnusedSpectra2;
                            } else {
                                lAverageInt2 = 1.0;
                            }

                            double lRest = Math.log( lAverageInt1 / lAverageInt2)/ Math.log(2);
                            if(Double.isNaN(lRest) || Double.isInfinite(lRest)){
                                //System.out.println(lRest);
                            } else {
                                //lCoRibarStat.addValue(lRest);
                            }
                        } else {
                            //lCoRibarStat.addValue(Math.log((lP1.getSI(true)/lP1.countNumberOfSpectra())/(lP2.getSI(true)/lP2.countNumberOfSpectra())) / Math.log(2));
                        }
                        */

                        double lExtraRibarValue = Math.log((lP1.getSI(true)/lP1.countNumberOfSpectra())/(lP2.getSI(true)/lP2.countNumberOfSpectra())) / Math.log(2);
                        lCoRibarStat.addValue(lExtraRibarValue);
                        //lCoRibarStat.addValue(Math.log((lP1.getSI(true)/(lP1.countNumberOfSpectra() - lCounter1))/(lP2.getSI(true)/(lP2.countNumberOfSpectra() - lCounter2))) / Math.log(2));
                        //lCoRibarStat.addValue(Math.log((lP1.getSI(true))/(lP2.getSI(true))) / Math.log(2));

                        double lNotUsedValues = (lIntSum1NotUsedInRibar / lCount1NotUsedInRibar) / (lIntSum2NotUsedInRibar / lCount2NotUsedInRibar);


                        double lRibar = 0.0;
                        double lCoRibar = 0.0;

                        lRibar = lRibarStat.getMean();
                        lRibar =  Math.pow(2,lRibar);
                        lCoRibar = lCoRibarStat.getMean();
                        lCoRibar =  Math.pow(2,lCoRibar);






                        lCounter = lCounter + 1;
                        double lEmPAIRatio = lProteins1.get(i).getEmPAI(true) / lProteins2.get(j).getEmPAI(true);
                        double lSINRatioMatched = lProteins1.get(i).getSIN(false) / lProteins2.get(j).getSIN(false);
                        double lNSAFRatio = lProteins1.get(i).getNSAF() / lProteins2.get(j).getNSAF();
                        double lSC = lProteins1.get(i).getSSSC() / lProteins2.get(j).getSSSC();

                        out.write(lProteins1.get(i).getAccession()+ "," + lStats1.size() + "," + lProteins1.get(i).countNumberOfSpectra() + ","  + lProteins1.get(i).countNumberOfUniqueParentIons()+ "," + lProteins2.get(j).countNumberOfSpectra() + ","  + lProteins2.get(j).countNumberOfUniqueParentIons() + ","  + lEmPAIRatio + "," + lSINRatioMatched + ","  + lNSAFRatio + "," + lSC + "," + lRibar + "," + lCoRibar + "," + lUPSstatus + "," + lRevString + "," + lExtraRibarValue + "," + lNotUsedValues + "\n");

                    }
                }
            }
            System.out.println(lCounter);
            out.flush();
            out.close();
        } catch(Exception e){
            e.printStackTrace();
        }





    }

}