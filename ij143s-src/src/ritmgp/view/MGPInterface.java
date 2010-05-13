package ritmgp.view;

import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.WindowManager;
import ij.gui.Plot;
import java.awt.Color;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;
import javax.swing.JFileChooser;
import ritmgp.model.MGPAlgorithm;
import ritmgp.model.MGPResult;
import ritmgp.util.MathUtil;

/**
 *
 * @author Mike Neurohr
 */
public class MGPInterface extends javax.swing.JFrame {

    private String[] imageTitles;
    private MGPResult result;

    public MGPInterface(String[] imageTitles) {
        initComponents();
        this.imageTitles = imageTitles;
        initComboBoxes();
        initTextFields();
    }

    public MGPInterface(List<String> imageTitles) {
        initComponents();
        this.imageTitles = imageTitles.toArray(new String[0]);
        initComboBoxes();
        initTextFields();
    }

    /** This method is called from within the constructor to
     * initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is
     * always regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        jScrollPane1 = new javax.swing.JScrollPane();
        resultsTable = new javax.swing.JTable();
        saveVecButton = new javax.swing.JButton();
        jLabel1 = new javax.swing.JLabel();
        jLabel2 = new javax.swing.JLabel();
        jLabel3 = new javax.swing.JLabel();
        jLabel4 = new javax.swing.JLabel();
        jLabel5 = new javax.swing.JLabel();
        jLabel6 = new javax.swing.JLabel();
        baseCheckBox = new javax.swing.JCheckBox();
        jLabel7 = new javax.swing.JLabel();
        leftTextField = new javax.swing.JTextField();
        jLabel8 = new javax.swing.JLabel();
        rightTextField = new javax.swing.JTextField();
        reflTextField = new javax.swing.JTextField();
        areaTextField = new javax.swing.JTextField();
        fStopTextField = new javax.swing.JTextField();
        fovTextField = new javax.swing.JTextField();
        dComboBox = new javax.swing.JComboBox();
        bComboBox = new javax.swing.JComboBox();
        jLabel9 = new javax.swing.JLabel();
        analyzeButton = new javax.swing.JButton();
        saveBRDFButton = new javax.swing.JButton();
        refreshButton = new javax.swing.JButton();
        jLabel10 = new javax.swing.JLabel();
        jLabel11 = new javax.swing.JLabel();
        jLabel12 = new javax.swing.JLabel();

        setDefaultCloseOperation(javax.swing.WindowConstants.EXIT_ON_CLOSE);
        setTitle("RITMGP Plugin");

        resultsTable.setModel(new javax.swing.table.DefaultTableModel(
            new Object [][] {
                {null, null, null, null, null, null, null, null, null, null}
            },
            new String [] {
                "area", "whalf", "w10", "h", "rho", "granularity", "aperture dim.", "sigma", "skewness", "kurtosis"
            }
        ));
        resultsTable.setEnabled(false);
        jScrollPane1.setViewportView(resultsTable);

        saveVecButton.setText("Save result vector");
        saveVecButton.setToolTipText("Save the above table");
        saveVecButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                saveVecButtonActionPerformed(evt);
            }
        });

        jLabel1.setText("Bright Image:");

        jLabel2.setText("Dark Image:");

        jLabel3.setText("Reference Reflectance:");

        jLabel4.setText("Reference Area:");

        jLabel5.setText("F-stop Measurement:");

        jLabel6.setText("Field of View:");

        baseCheckBox.setText("Auto Baseline");
        baseCheckBox.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                baseCheckBoxActionPerformed(evt);
            }
        });

        jLabel7.setText("Minimum angle:");

        jLabel8.setText("Maximum angle:");

        dComboBox.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "Item 1", "Item 2", "Item 3", "Item 4" }));

        bComboBox.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "Item 1", "Item 2", "Item 3", "Item 4" }));

        jLabel9.setText("Enter 1 if measuring a reference sample");

        analyzeButton.setText("Analyze");
        analyzeButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                analyzeButtonActionPerformed(evt);
            }
        });

        saveBRDFButton.setText("Save BRDF");
        saveBRDFButton.setToolTipText("Save the final, baseline-corrected BRDF");
        saveBRDFButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                saveBRDFButtonActionPerformed(evt);
            }
        });

        refreshButton.setText("Refresh images");
        refreshButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                refreshButtonActionPerformed(evt);
            }
        });

        jLabel10.setText("mm");

        jLabel11.setText("deg");

        jLabel12.setText("deg");

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(getContentPane());
        getContentPane().setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                    .addComponent(jScrollPane1, javax.swing.GroupLayout.DEFAULT_SIZE, 609, Short.MAX_VALUE)
                    .addGroup(layout.createSequentialGroup()
                        .addComponent(saveBRDFButton)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(saveVecButton))
                    .addGroup(javax.swing.GroupLayout.Alignment.LEADING, layout.createSequentialGroup()
                        .addComponent(jLabel7)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(leftTextField, javax.swing.GroupLayout.PREFERRED_SIZE, 42, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(jLabel11)
                        .addGap(18, 18, 18)
                        .addComponent(jLabel8)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(rightTextField, javax.swing.GroupLayout.PREFERRED_SIZE, 43, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(jLabel12))
                    .addComponent(analyzeButton, javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(javax.swing.GroupLayout.Alignment.LEADING, layout.createSequentialGroup()
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(jLabel3)
                            .addComponent(jLabel2)
                            .addComponent(jLabel1)
                            .addComponent(jLabel5)
                            .addComponent(jLabel6)
                            .addComponent(jLabel4))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                            .addComponent(reflTextField)
                            .addComponent(dComboBox, 0, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(bComboBox, 0, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(areaTextField)
                            .addComponent(fStopTextField)
                            .addComponent(fovTextField))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(jLabel9)
                            .addComponent(refreshButton)
                            .addComponent(jLabel10)))
                    .addComponent(baseCheckBox, javax.swing.GroupLayout.Alignment.LEADING))
                .addContainerGap())
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jLabel1)
                    .addComponent(bComboBox, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jLabel2)
                    .addComponent(dComboBox, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(refreshButton))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jLabel3)
                    .addComponent(reflTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(areaTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(jLabel4)
                    .addComponent(jLabel9))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jLabel5)
                    .addComponent(fStopTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jLabel6)
                    .addComponent(fovTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(jLabel10))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(baseCheckBox)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jLabel7)
                    .addComponent(leftTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(jLabel11)
                    .addComponent(rightTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(jLabel8)
                    .addComponent(jLabel12))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(analyzeButton)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addComponent(jScrollPane1, javax.swing.GroupLayout.PREFERRED_SIZE, 40, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addGap(8, 8, 8)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(saveVecButton)
                    .addComponent(saveBRDFButton))
                .addContainerGap())
        );

        pack();
    }// </editor-fold>//GEN-END:initComponents

    private void analyzeButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_analyzeButtonActionPerformed
        if (imageTitles.length == 0 || bComboBox.getSelectedItem() == null
                || dComboBox.getSelectedItem() == null) {
            IJ.showMessage("Missing at least one image selection.");
            return;
        }
        if(bComboBox.getSelectedItem().toString().equals(
                dComboBox.getSelectedItem().toString())){
            IJ.showMessage("Bright and dark images are the same image");
            return;
        }
        try {
            ImagePlus bright = WindowManager.getImage(bComboBox.getSelectedItem().toString());
            ImagePlus dark = WindowManager.getImage(dComboBox.getSelectedItem().toString());
            double refRefl = Double.parseDouble(reflTextField.getText());
            double refArea = Double.parseDouble(areaTextField.getText());
            double fStop = Double.parseDouble(fStopTextField.getText());
            double fov = Double.parseDouble(fovTextField.getText());
            boolean baseline = baseCheckBox.isSelected();
            double leftAngle = baseline ? 0 : Double.parseDouble(leftTextField.getText());
            double rightAngle = baseline ? 0 : Double.parseDouble(rightTextField.getText());
            MGPAlgorithm algo = new MGPAlgorithm(bright, dark, refArea, refRefl,
                    fov, fStop, baseline, leftAngle, rightAngle);
            result = algo.measureGloss();
            double[] xs = result.getAlpha();
            Plot BRDF = new Plot("BRDF", "alpha", "intensity", xs, result.getGraphBRDF());
            BRDF.setLimits(MathUtil.min(xs), MathUtil.max(xs),
                    MathUtil.min(result.getGraphBRDF()) - MathUtil.max(result.getGraphBRDF()) / 5,
                    MathUtil.max(result.getGraphBRDF()) + MathUtil.max(result.getGraphBRDF()) / 5);
            double baselineX1 = result.getLeftAngle();
            double baselineY1 = MathUtil.linterp(xs, result.getGraphBRDF(), baselineX1);
            double baselineX2 = result.getRightAngle();
            double baselineY2 = MathUtil.linterp(xs, result.getGraphBRDF(), baselineX2);
            BRDF.draw();
            BRDF.setColor(Color.red);
            BRDF.drawLine(baselineX1, baselineY1, baselineX2, baselineY2);
            BRDF.show();
            resultsTable.setModel(new javax.swing.table.DefaultTableModel(
                    new Object[][]{
                        result.getVector()
                    },
                    new String[]{
                        "area", "whalf", "w10", "h", "rho", "granularity", "aperture dim.", "sigma", "skewness", "kurtosis"
                    }));
            Prefs.set("RITMGP.refl", reflTextField.getText());
            Prefs.set("RITMGP.area", areaTextField.getText());
            Prefs.set("RITMGP.exp", fStopTextField.getText());
            Prefs.set("RITMGP.fov", fovTextField.getText());
            Prefs.set("RITMGP.left", leftTextField.isEditable()
                    ? leftTextField.getText() : "");
            Prefs.set("RITMGP.right", rightTextField.isEditable()
                    ? rightTextField.getText() : "");
            Prefs.set("RITMGP.base", baseline);
            Prefs.savePreferences();
        } catch (NumberFormatException e) {
            IJ.showMessage("Please enter a value in all editable fields.");
        }
    }//GEN-LAST:event_analyzeButtonActionPerformed

    private void baseCheckBoxActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_baseCheckBoxActionPerformed
        leftTextField.setEditable(!baseCheckBox.isSelected());
        rightTextField.setEditable(!baseCheckBox.isSelected());
        leftTextField.setText("");
        rightTextField.setText("");
    }//GEN-LAST:event_baseCheckBoxActionPerformed

    private void saveVecButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_saveVecButtonActionPerformed
        if (result == null) {
            IJ.showMessage("No results yet!");
            return;
        }
        JFileChooser fileChooser = new JFileChooser();
        fileChooser.setMultiSelectionEnabled(false);
        if (fileChooser.showSaveDialog(this) == JFileChooser.APPROVE_OPTION) {
            PrintWriter vecFile = null;
            try {
                vecFile = new PrintWriter(new BufferedWriter(new FileWriter(
                        fileChooser.getSelectedFile())));
                String[] titles = {"area", "whalf", "w10", "h", "rho",
                    "granularity", "aperture dim.", "sigma", "skewness", "kurtosis"};
                for (int i = 0; i < result.getVector().length; i++) {
                    vecFile.println(titles[i] + "\t" + result.getVector()[i]);
                }
            } catch (IOException ex) {
                IJ.showMessage("Error writing to file, "
                        + fileChooser.getSelectedFile().getName() + ": "
                        + ex.getMessage());
            } finally {
                if (vecFile != null) {
                    vecFile.close();
                }
            }

        }

    }//GEN-LAST:event_saveVecButtonActionPerformed

    private void saveBRDFButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_saveBRDFButtonActionPerformed
        if (result == null) {
            IJ.showMessage("No results yet!");
            return;
        }
        JFileChooser fileChooser = new JFileChooser();
        fileChooser.setMultiSelectionEnabled(false);
        if (fileChooser.showSaveDialog(this) == JFileChooser.APPROVE_OPTION) {
            PrintWriter brdfFile = null;
            try {
                brdfFile = new PrintWriter(new BufferedWriter(new FileWriter(
                        fileChooser.getSelectedFile())));
                for (int i = 0; i < result.getGraphBRDF().length; i++) {
                    brdfFile.println(result.getAlpha()[i] + "\t" + result.getFinalBRDF()[i]);
                }
            } catch (IOException ex) {
                IJ.showMessage("Error writing to file, "
                        + fileChooser.getSelectedFile().getName() + ": "
                        + ex.getMessage());
            } finally {
                if (brdfFile != null) {
                    brdfFile.close();
                }
            }

        }
    }//GEN-LAST:event_saveBRDFButtonActionPerformed

    private void refreshButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_refreshButtonActionPerformed
        refreshImages();
    }//GEN-LAST:event_refreshButtonActionPerformed
    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JButton analyzeButton;
    private javax.swing.JTextField areaTextField;
    private javax.swing.JComboBox bComboBox;
    private javax.swing.JCheckBox baseCheckBox;
    private javax.swing.JComboBox dComboBox;
    private javax.swing.JTextField fStopTextField;
    private javax.swing.JTextField fovTextField;
    private javax.swing.JLabel jLabel1;
    private javax.swing.JLabel jLabel10;
    private javax.swing.JLabel jLabel11;
    private javax.swing.JLabel jLabel12;
    private javax.swing.JLabel jLabel2;
    private javax.swing.JLabel jLabel3;
    private javax.swing.JLabel jLabel4;
    private javax.swing.JLabel jLabel5;
    private javax.swing.JLabel jLabel6;
    private javax.swing.JLabel jLabel7;
    private javax.swing.JLabel jLabel8;
    private javax.swing.JLabel jLabel9;
    private javax.swing.JScrollPane jScrollPane1;
    private javax.swing.JTextField leftTextField;
    private javax.swing.JTextField reflTextField;
    private javax.swing.JButton refreshButton;
    private javax.swing.JTable resultsTable;
    private javax.swing.JTextField rightTextField;
    private javax.swing.JButton saveBRDFButton;
    private javax.swing.JButton saveVecButton;
    // End of variables declaration//GEN-END:variables

    private void initComboBoxes() {
        bComboBox.setModel(new javax.swing.DefaultComboBoxModel(imageTitles));
        dComboBox.setModel(new javax.swing.DefaultComboBoxModel(imageTitles));
    }

    private void refreshImages() {
        int[] wList = WindowManager.getIDList(); //ij.WindowManager
        if (wList == null) {
            imageTitles = new String[0];
            initComboBoxes();
            return;
        }
        ArrayList<String> titles = new ArrayList<String>(wList.length);
        int i = 0;
        for (int imgId : wList) {
            ImagePlus imp = WindowManager.getImage(imgId);
            if (imp != null) {
                titles.add(imp.getTitle());
            }
        }
        imageTitles = titles.toArray(new String[0]);
        initComboBoxes();
    }

    private void initTextFields() {
        reflTextField.setText(Prefs.get("RITMGP.refl", ""));
        areaTextField.setText(Prefs.get("RITMGP.area", ""));
        fStopTextField.setText(Prefs.get("RITMGP.exp", ""));
        fovTextField.setText(Prefs.get("RITMGP.fov", ""));
        baseCheckBox.setSelected(Prefs.getBoolean("RITMGP.base", false));
        leftTextField.setText(Prefs.get("RITMGP.left", ""));
        rightTextField.setText(Prefs.get("RITMGP.right", ""));
        baseCheckBox.setSelected(Prefs.getBoolean("RITMGP.base", false));
        if (baseCheckBox.isSelected()) {
            leftTextField.setEditable(false);
            rightTextField.setEditable(false);
            leftTextField.setText("");
            rightTextField.setText("");
        }
    }
}
