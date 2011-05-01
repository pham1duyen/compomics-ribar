package com.compomics.ribar;

import com.compomics.mascotdatfile.util.interfaces.MascotDatfileInf;
import com.compomics.mascotdatfile.util.interfaces.QueryToPeptideMapInf;
import com.compomics.mascotdatfile.util.mascot.*;
import com.compomics.mascotdatfile.util.mascot.fragmentions.FragmentIonImpl;
import com.compomics.mslims.db.accessors.Identification;
import com.compomics.mslims.db.accessors.IdentificationTableAccessor;
import com.compomics.mslims.util.mascot.MascotIdentifiedSpectrum;
import com.compomics.mslims.util.mascot.MascotIsoforms;
import com.sun.java.swing.plaf.nimbus.NimbusLookAndFeel;

import javax.swing.*;
import javax.swing.filechooser.FileFilter;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.*;
import java.math.BigDecimal;
import java.sql.SQLException;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: Niklaas
 * Date: 29/04/11
 * Time: 12:50
 * To change this template use File | Settings | File Templates.
 */
public class RibarGui extends JFrame {
    private JTextField txtOutFile;
    private JButton iOpenButton;
    private JButton iOpenButton3;
    private JButton iOpenButton2;
    private JLabel lblSet1Count;
    private JButton iRunButton;
    private JProgressBar iProgressBar1;
    private JTextField txtDatFileValue;
    private JLabel lblSet2Count;
    private JPanel contentPane;
    private JButton iCompOmics;
    private Double iDatfileValue = 0.05;
    public Vector<File> iSetOneFiles = new Vector<File>();
    public Vector<File> iSetTwoFiles = new Vector<File>();
    public Vector<MascotIdentifiedSpectrumExtension> iSetOneIdentifications = new Vector<MascotIdentifiedSpectrumExtension>();
    public Vector<MascotIdentifiedSpectrumExtension> iSetTwoIdentifications = new Vector<MascotIdentifiedSpectrumExtension>();
    private File iOutFile;
    private MascotDatfile_Index iDatfile;


    public RibarGui() {
        this.setTitle("RIBAR and xRIBAR calculations");
        this.setIconImage(new ImageIcon(getClass().getResource("/compomics.png")).getImage());

        txtDatFileValue.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                String lMsMsError = txtDatFileValue.getText();
                try {
                    iDatfileValue = Double.valueOf(lMsMsError);
                } catch (NumberFormatException es) {
                    txtDatFileValue.setText(String.valueOf(iDatfileValue));
                    JOptionPane.showMessageDialog(getFrame(), iDatfileValue + " is not a valid value!", "Number error", JOptionPane.WARNING_MESSAGE);
                }

            }
        });
        iOpenButton3.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                //open file chooser
                JFileChooser fc = new JFileChooser();
                fc.setMultiSelectionEnabled(true);
                //create the file filter to choose
                FileFilter lFilter = new DatFileFilter();
                fc.setFileFilter(lFilter);
                fc.showOpenDialog(getFrame());
                File[] lFiles = fc.getSelectedFiles();
                iSetOneFiles.removeAllElements();
                for (int i = 0; i < lFiles.length; i++) {
                    iSetOneFiles.add(lFiles[i]);
                }
                lblSet1Count.setText(iSetOneFiles.size() + " .dat files added to set one");
                testEverythingSet();
            }
        });
        iOpenButton2.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                //open file chooser
                JFileChooser fc = new JFileChooser();
                fc.setMultiSelectionEnabled(true);
                //create the file filter to choose
                FileFilter lFilter = new DatFileFilter();
                fc.setFileFilter(lFilter);
                fc.showOpenDialog(getFrame());
                File[] lFiles = fc.getSelectedFiles();
                iSetTwoFiles.removeAllElements();
                for (int i = 0; i < lFiles.length; i++) {
                    iSetTwoFiles.add(lFiles[i]);
                }
                lblSet2Count.setText(iSetTwoFiles.size() + " .dat files added to set two");
                testEverythingSet();
            }
        });
        iOpenButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                //open file chooser
                JFileChooser fc = new JFileChooser();
                fc.setMultiSelectionEnabled(false);
                FileFilter lFilter = new CsvFileFilter();
                fc.setFileFilter(lFilter);
                fc.showOpenDialog(getFrame());
                iOutFile = fc.getSelectedFile();
                if (!iOutFile.getAbsolutePath().endsWith(".csv")) {
                    iOutFile = new File(iOutFile.getAbsolutePath() + ".csv");
                }
                if (iOutFile != null) {
                    txtOutFile.setText(iOutFile.getAbsolutePath());
                }
                testEverythingSet();
            }
        });
        iRunButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                run();
            }
        });
        //add a closing window listener
        addWindowListener(new WindowAdapter() {
            public void windowClosing(WindowEvent evt) {
                System.exit(0);
            }
        }
        );
        iCompOmics.setIcon(new ImageIcon(getClass().getResource("/compomics.png")));
        iCompOmics.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                System.out.println("http://www.compomics.com");
                showInBrowser("http://www.compomics.com");
            }
        });

        iProgressBar1.setVisible(false);
        this.setContentPane(contentPane);
        this.setLocationRelativeTo(null);
        this.setSize(540, 350);
        this.setVisible(true);
        Dimension d = Toolkit.getDefaultToolkit().getScreenSize();
        int x = (d.width - getSize().width) / 10;
        int y = (d.height - getSize().height) / 10;
        this.setLocation(x, y);
    }

    public void setClickable(boolean lClicksAllowed) {
        iOpenButton.setEnabled(lClicksAllowed);
        iOpenButton2.setEnabled(lClicksAllowed);
        iOpenButton3.setEnabled(lClicksAllowed);
        iRunButton.setEnabled(lClicksAllowed);
        txtDatFileValue.setEnabled(lClicksAllowed);
    }

    /**
     * This method opens the default browser on a given webpage
     *
     * @param url String with the url
     * @return boolean False if an error occured
     */
    private boolean showInBrowser(String url) {

        String os = System.getProperty("os.name").toLowerCase();
        Runtime rt = Runtime.getRuntime();
        try {
            if (os.indexOf("win") >= 0) {
                String[] cmd = new String[4];
                cmd[0] = "cmd.exe";
                cmd[1] = "/C";
                cmd[2] = "start";
                cmd[3] = url;
                rt.exec(cmd);
            } else if (os.indexOf("mac") >= 0) {
                rt.exec("open " + url);
            } else {
                //prioritized 'guess' of users' preference
                String[] browsers = {"epiphany", "firefox", "mozilla", "konqueror",
                        "netscape", "opera", "links", "lynx"};

                StringBuffer cmd = new StringBuffer();
                for (int i = 0; i < browsers.length; i++)
                    cmd.append((i == 0 ? "" : " || ") + browsers[i] + " \"" + url + "\" ");

                rt.exec(new String[]{"sh", "-c", cmd.toString()});
                //rt.exec("firefox http://www.google.com");
                //System.out.println(cmd.toString());

            }
        } catch (IOException e) {
            JOptionPane.showMessageDialog(new JFrame(), "\n\n The system failed to invoke your default web browser while attempting to access: \n\n " + url + "\n\n", "Browser Error", JOptionPane.WARNING_MESSAGE);
            return false;
        }
        return true;
    }

    public void run() {
        com.compomics.util.sun.SwingWorker lRunner = new com.compomics.util.sun.SwingWorker() {
            public Boolean construct() {
                setClickable(false);
                iProgressBar1.setStringPainted(true);
                iProgressBar1.setVisible(true);
                iProgressBar1.setMaximum(iSetOneFiles.size() + iSetTwoFiles.size() + 1);
                //get the identification from set 1
                for (int i = 0; i < iSetOneFiles.size(); i++) {
                    iProgressBar1.setValue(i + 1);
                    iProgressBar1.setString("Parsing datfile " + (i + 1) + " of  " + iSetOneFiles.size() + " in set one.");

                    try {
                        iDatfile = new MascotDatfile_Index(new BufferedReader(new FileReader(iSetOneFiles.get(i))), iSetOneFiles.get(i).getName());
                        Vector<MascotIdentifiedSpectrumExtension> lIdentifications = extractIDs(iDatfile);
                        for (int k = 0; k < lIdentifications.size(); k++) {
                            iSetOneIdentifications.add(lIdentifications.get(k));
                        }

                    } catch (FileNotFoundException e) {
                        e.printStackTrace();
                    }
                }


                for (int i = 0; i < iSetTwoFiles.size(); i++) {
                    iProgressBar1.setValue(iSetOneFiles.size() + i + 1);
                    iProgressBar1.setString("Parsing datfile " + (i + 1) + " of  " + iSetTwoFiles.size() + " in set two.");

                    try {
                        iDatfile = new MascotDatfile_Index(new BufferedReader(new FileReader(iSetTwoFiles.get(i))), iSetTwoFiles.get(i).getName());
                        Vector<MascotIdentifiedSpectrumExtension> lIdentifications = extractIDs(iDatfile);
                        for (int k = 0; k < lIdentifications.size(); k++) {
                            iSetTwoIdentifications.add(lIdentifications.get(k));
                        }

                    } catch (FileNotFoundException e) {
                        e.printStackTrace();
                    }
                }

                iProgressBar1.setIndeterminate(true);
                iProgressBar1.setString("Collecting all information!");
                Project lProject1 = new Project(iSetOneIdentifications);
                Project lProject2 = new Project(iSetTwoIdentifications);

                iProgressBar1.setString("Writing csv file!");
                new SimpleComparer(lProject1, lProject2, iOutFile.getAbsolutePath());
                iProgressBar1.setIndeterminate(false);
                iSetOneFiles.removeAllElements();
                iSetTwoFiles.removeAllElements();
                iOutFile = null;
                testEverythingSet();

                JOptionPane.showMessageDialog(getFrame(), "All data was saved to the csv file.\nCalculated RIBAR and xRIBAR (set 1 / set 2) protein ratios.", "RIBAR and xRIBAR completed", JOptionPane.INFORMATION_MESSAGE);
                setClickable(true);

                return true;
            }

            public void finished() {
                iProgressBar1.setVisible(false);
            }

        };
        lRunner.start();

    }


    private Vector<MascotIdentifiedSpectrumExtension> extractIDs(MascotDatfileInf aMDF) throws IllegalArgumentException {
        // Vector that will contain the MascotIdentifiedSpectrum instances.
        Vector<MascotIdentifiedSpectrumExtension> result = new Vector<MascotIdentifiedSpectrumExtension>();

        // Get the generic parameters for the search,
        // Extract the db filename and the Mascot version.
        Header header = aMDF.getHeaderSection();
        String version = header.getVersion();
        String dbfilename = header.getRelease();
        Parameters parameters = aMDF.getParametersSection();
        String searchTitle = parameters.getCom();
        if (searchTitle == null) {
            searchTitle = "!No title specified";
        } else {
            int location = searchTitle.indexOf("|");
            if (location >= 0) {
                searchTitle = searchTitle.substring(0, location).trim();
            }
        }
        String inputfile = parameters.getFile();
        String dbName = parameters.getDatabase();
        ProteinMap proteinMap = aMDF.getProteinMap();
        Masses masses = aMDF.getMasses();

        // Rank of the hit (only highest ranking hits
        // (i.e.: rank = 1)) are considered,
        int rank = 1;

        System.out.println("Get querytopepmap");
        QueryToPeptideMapInf queryToPepMap = aMDF.getQueryToPeptideMap();
        System.out.println("Start iterator");
        System.out.println("Get query list");
        // Get all the queries...
        int lQueryCounter = 0;
        for (int i = 1; i <= aMDF.getNumberOfQueries(); i++) {

            if (lQueryCounter % 1000 == 0) {
                System.out.println(lQueryCounter + " / " + aMDF.getNumberOfQueries());
            }
            // Get the query.
            Query query = aMDF.getQuery(i);
            lQueryCounter++;
            // Get the first ranking peptide hit, if any.
            PeptideHit ph = queryToPepMap.getPeptideHitOfOneQuery(query.getQueryNumber(), rank);
            if (ph != null && ph.scoresAboveIdentityThreshold(iDatfileValue)) {
                // We have a peptide hit for this query that scores equal
                // to or above the threshold. Parse it and create a
                // MascotIdentifiedSpectrum.
                Peak[] lPeaks = query.getPeakList();
                double lSummedPeaks = 0.0;
                for (int j = 0; j < lPeaks.length; j++) {
                    lSummedPeaks = lSummedPeaks + lPeaks[j].getIntensity();
                }
                MascotIdentifiedSpectrumExtension mis = new MascotIdentifiedSpectrumExtension(lSummedPeaks);
                // Generic stuff, already parsed in advance.
                mis.setDBFilename(dbfilename);
                mis.setMascotVersion(version);
                mis.setSearchTitle(searchTitle);
                mis.setOriginal_file(inputfile);
                mis.setDBName(dbName);
                mis.setQueryNr(lQueryCounter);

                // Query title.


                //if it is a multifile get the scans for the query
                mis.setFile(query.getTitle());

                // Additional query info.
                if (mis.getFile() == null && aMDF.getNumberOfQueries() == 1) {
                    // In this case, a single query was performed using a file that did not contain
                    // 'merge' (regardless of case). This is indicative of a search with a single spectrum.
                    // Therefore we just keep the name of the spectrum as reported by Mascot.
                    mis.setFile(inputfile);
                } else if (mis.getFile() == null && aMDF.getNumberOfQueries() > 1) {
                    // Mergefile.
                    // We omit the filename (set it to '*').
                    mis.setFile("*");
                }

                // Query m/z and charge.
                double mz = new BigDecimal(query.getPrecursorMZ()).setScale(4, BigDecimal.ROUND_HALF_UP).doubleValue();
                String chargeString = query.getChargeString();
                boolean isNegative = false;
                int chargeLoc = chargeString.indexOf('+');
                if (chargeLoc < 0) {
                    chargeLoc = chargeString.indexOf('-');
                    isNegative = true;
                }
                chargeString = chargeString.substring(0, chargeLoc);
                int charge = Integer.parseInt(chargeString);
                if (isNegative) {
                    charge = -charge;
                }
                mis.setChargeState(charge);
                mis.setPrecursorMZ(mz);

                // PeptideHit stuff.
                // Thresholds and rank.
                mis.setHomologyTreshold((int) ph.getHomologyThreshold());
                mis.setIdentityTreshold((int) ph.calculateIdentityThreshold(iDatfileValue));
                mis.setRank(rank);

                mis.setTheoreticalMass(ph.getPeptideMr());
                mis.setMeasuredMass(ph.getPeptideMr() + ph.getDeltaMass());
                mis.setSequence(ph.getSequence());
                String lModifiedSequence = ph.getModifiedSequence();
                // If a modified sequence contain's a '#' character, this means the modification was not included in the modificationConversion.txt file.
                // Throw an error since we don't want to have multiple names for identical modifications.
                if (lModifiedSequence.indexOf('#') != -1) {
                    //throw new IllegalArgumentException("\n\nModificationConversion.txt does not contain enough information to parse the following identification:\n\t" + lModifiedSequence + "\nPlease add the modification into modificationcoverions.txt. ");
                }
                mis.setModifiedSequence(lModifiedSequence);
                mis.setScore(ph.getIonsScore());

                // Protein stuff.
                MascotIsoforms mifs = new MascotIsoforms();
                Iterator iter2 = ph.getProteinHits().iterator();
                while (iter2.hasNext()) {
                    ProteinHit protein = (ProteinHit) iter2.next();
                    // Hold the original accession to access
                    String originalAccession = protein.getAccession();
                    String trimmedAccession = originalAccession;
                    int startLoc = trimmedAccession.indexOf('(');
                    int endLoc = trimmedAccession.indexOf(')');
                    int tempStart = -1;
                    int tempEnd = -1;
                    if ((startLoc >= 0) && (endLoc >= 0)) {
                        String tempLocalization = trimmedAccession.substring(startLoc + 1, endLoc);
                        StringTokenizer lst = new StringTokenizer(tempLocalization, "-");
                        try {
                            tempStart = Integer.parseInt(lst.nextToken().trim());
                            tempEnd = Integer.parseInt(lst.nextToken().trim());
                            trimmedAccession = trimmedAccession.substring(0, startLoc).trim();
                        } catch (Exception e) {
                            // Do nothing.
                            // It's probably just not a location String.
                        }
                    }
                    // If no start and end location found, take those from the
                    // protein information supplied by Mascot.
                    if (tempStart < 0) {
                        tempStart = protein.getStart();
                        tempEnd = protein.getStop();
                    }
                    mifs.addIsoform(trimmedAccession, proteinMap.getProteinDescription(originalAccession), tempStart, tempEnd);
                }
                mis.setIsoforms(mifs);
                // Add the ion coverage String.
                PeptideHitAnnotation pha = ph.getPeptideHitAnnotation(masses, parameters, query.getPrecursorMZ(), query.getChargeString());
                pha.getMatchedBYions(query.getPeakList());
                // Calling this method will initialize the ion importance as determined by Mascot.
                Collection fragmentions = pha.getFusedMatchedIons(query.getPeakList(), ph.getPeaksUsedFromIons1(), query.getMaxIntensity(), 0.10);
                mis.setFragmentIons(fragmentions);
                double fragmentError = Double.parseDouble(parameters.getITOL());
                String fragmentErrorUnit = parameters.getITOLU();
                if (fragmentErrorUnit.trim().toLowerCase().equals("ppm")) {
                    fragmentError = query.getPrecursorMZ() * fragmentError * 1e-6;
                }
                mis.setFragmentMassError(fragmentError);
                // Add mis to vector.
                result.add(mis);
            }
        }
        return result;
    }

    public void testEverythingSet() {
        if (iSetOneFiles.size() > 0 && iSetTwoFiles.size() > 0 && iOutFile != null) {
            iRunButton.setEnabled(true);
        } else {
            iRunButton.setEnabled(false);
        }
    }

    public JFrame getFrame() {
        return this;
    }


    public static void main(String[] args) {
        try {
            UIManager.setLookAndFeel(new NimbusLookAndFeel());
        } catch (UnsupportedLookAndFeelException e) {
            // ignore exception
        }
        new RibarGui();
    }

    {
// GUI initializer generated by IntelliJ IDEA GUI Designer
// >>> IMPORTANT!! <<<
// DO NOT EDIT OR ADD ANY CODE HERE!
        $$$setupUI$$$();
    }

    /**
     * Method generated by IntelliJ IDEA GUI Designer
     * >>> IMPORTANT!! <<<
     * DO NOT edit this method OR call it in your code!
     *
     * @noinspection ALL
     */
    private void $$$setupUI$$$() {
        contentPane = new JPanel();
        contentPane.setLayout(new GridBagLayout());
        final JLabel label1 = new JLabel();
        label1.setText("Csv output file: ");
        GridBagConstraints gbc;
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 3;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.insets = new Insets(5, 5, 5, 5);
        contentPane.add(label1, gbc);
        txtOutFile = new JTextField();
        txtOutFile.setEditable(false);
        gbc = new GridBagConstraints();
        gbc.gridx = 1;
        gbc.gridy = 3;
        gbc.weightx = 1.0;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(5, 5, 5, 5);
        contentPane.add(txtOutFile, gbc);
        iOpenButton = new JButton();
        iOpenButton.setText("Open");
        gbc = new GridBagConstraints();
        gbc.gridx = 2;
        gbc.gridy = 3;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(5, 5, 5, 5);
        contentPane.add(iOpenButton, gbc);
        final JLabel label2 = new JLabel();
        label2.setText("Set 1 Mascot .dat files");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 0;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.insets = new Insets(5, 5, 5, 5);
        contentPane.add(label2, gbc);
        final JLabel label3 = new JLabel();
        label3.setText("Set 2 Mascot .dat files");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 1;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.insets = new Insets(5, 5, 5, 5);
        contentPane.add(label3, gbc);
        final JLabel label4 = new JLabel();
        label4.setText("Mascot .dat file alpha-value");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 2;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.insets = new Insets(5, 5, 5, 5);
        contentPane.add(label4, gbc);
        iOpenButton3 = new JButton();
        iOpenButton3.setText("Open");
        gbc = new GridBagConstraints();
        gbc.gridx = 2;
        gbc.gridy = 0;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(5, 5, 5, 5);
        contentPane.add(iOpenButton3, gbc);
        iOpenButton2 = new JButton();
        iOpenButton2.setText("Open");
        gbc = new GridBagConstraints();
        gbc.gridx = 2;
        gbc.gridy = 1;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(5, 5, 5, 5);
        contentPane.add(iOpenButton2, gbc);
        lblSet1Count = new JLabel();
        lblSet1Count.setFont(new Font(lblSet1Count.getFont().getName(), Font.ITALIC, lblSet1Count.getFont().getSize()));
        lblSet1Count.setText(" No files added to set one");
        gbc = new GridBagConstraints();
        gbc.gridx = 1;
        gbc.gridy = 0;
        gbc.insets = new Insets(5, 5, 5, 5);
        contentPane.add(lblSet1Count, gbc);
        lblSet2Count = new JLabel();
        lblSet2Count.setFont(new Font(lblSet2Count.getFont().getName(), Font.ITALIC, lblSet2Count.getFont().getSize()));
        lblSet2Count.setText(" No files added to set two");
        gbc = new GridBagConstraints();
        gbc.gridx = 1;
        gbc.gridy = 1;
        gbc.insets = new Insets(5, 5, 5, 5);
        contentPane.add(lblSet2Count, gbc);
        final JSeparator separator1 = new JSeparator();
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 4;
        gbc.gridwidth = 3;
        gbc.fill = GridBagConstraints.BOTH;
        gbc.insets = new Insets(5, 5, 5, 5);
        contentPane.add(separator1, gbc);
        iRunButton = new JButton();
        iRunButton.setEnabled(false);
        iRunButton.setText("Run");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 5;
        gbc.gridwidth = 3;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(5, 5, 5, 5);
        contentPane.add(iRunButton, gbc);
        iProgressBar1 = new JProgressBar();
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 6;
        gbc.gridwidth = 3;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(5, 5, 5, 5);
        contentPane.add(iProgressBar1, gbc);
        txtDatFileValue = new JTextField();
        txtDatFileValue.setHorizontalAlignment(4);
        txtDatFileValue.setText("0.05");
        gbc = new GridBagConstraints();
        gbc.gridx = 2;
        gbc.gridy = 2;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(5, 5, 5, 5);
        contentPane.add(txtDatFileValue, gbc);
        iCompOmics = new JButton();
        iCompOmics.setBorderPainted(false);
        iCompOmics.setContentAreaFilled(false);
        iCompOmics.setFocusPainted(false);
        iCompOmics.setRolloverEnabled(false);
        iCompOmics.setText("");
        gbc = new GridBagConstraints();
        gbc.gridx = 2;
        gbc.gridy = 7;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        contentPane.add(iCompOmics, gbc);
    }

    /**
     * @noinspection ALL
     */
    public JComponent $$$getRootComponent$$$() {
        return contentPane;
    }


    /**
     * A .dat file filter
     */
    class DatFileFilter extends FileFilter {
        public boolean accept(File f) {
            return f.isDirectory() || f.getName().toLowerCase().endsWith(".dat");
        }

        public String getDescription() {
            return ".dat files";
        }
    }

    /**
     * A .csv file filter
     */
    class CsvFileFilter extends FileFilter {
        public boolean accept(File f) {
            return f.isDirectory() || f.getName().toLowerCase().endsWith(".csv");
        }

        public String getDescription() {
            return ".csv files";
        }
    }
}
