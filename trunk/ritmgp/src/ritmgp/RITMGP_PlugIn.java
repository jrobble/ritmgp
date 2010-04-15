package ritmgp;

import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.plugin.PlugIn;
import java.util.ArrayList;
import javax.swing.JFrame;
import ritmgp.view.MGPInterface;

/**
 *
 * @author Mike Neurohr
 */
public class RITMGP_PlugIn implements PlugIn{

    public void run(String string) {
       int[] wList = WindowManager.getIDList(); //ij.WindowManager
       if (wList==null || wList.length<2) {
           IJ.showMessage("Error", "There must be at least two windows open"); //ij.IJ
           return;
       }
       ArrayList<String> titles = new ArrayList<String>(wList.length);
       int i = 0;
       for (int imgId : wList) {
           ImagePlus imp = WindowManager.getImage(imgId);
           if (imp!=null)
               titles.add(imp.getTitle());
       }

       MGPInterface inter = new MGPInterface(titles);
       inter.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
       inter.setVisible(true);
    }

}
