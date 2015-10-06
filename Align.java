/*  Align:  main class

    Alignment indexes: left-most index (on the optical map or the forward strand
    of the reference maps), starting from 0, all inclusive.
    
    
    --
    OPTIMA v.f-1.3 -- 6 October 2015
    Copyright (C) Davide Verzotto and Niranjan Nagarajan

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License 2.1 as published by the Free Software Foundation.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License 2.1 along with this library; if not, write to the Free Software
    Foundation, Inc., 51 Franklin Street, 5th Floor, Boston, MA 02110-1301, USA.

 */
 
 
package OPTIMA;

import java.io.IOException;
import java.io.PrintWriter;
import java.io.FileWriter;

public class Align
{
  
  public static void main (String[] arg) throws IOException
  { if (arg.length < 5)
    { System.err.println("Use: java OPTIMA/Align.java OpticalMaps.maps OutputFileName InSilicoMaps.silico [pvalue|score] [allMaps|select] (OPTIONAL: firstIndex lastIndex)");
      System.exit(0);
    }

    long timeIn = System.currentTimeMillis();

    ExtendSeeds program1ExtendSeeds = new ExtendSeeds();


    program1ExtendSeeds.out = new PrintWriter(new FileWriter(arg[1] + ".txt"), true);
    program1ExtendSeeds.outOK = new PrintWriter(new FileWriter(arg[1] + ".ok"), true);
    program1ExtendSeeds.outNotOK = new PrintWriter(new FileWriter(arg[1] + ".notOk"), true);
    program1ExtendSeeds.outNotFound = new PrintWriter(new FileWriter(arg[1] + ".notFound"), true);
    program1ExtendSeeds.outDiscarded = new PrintWriter(new FileWriter(arg[1] + ".discarded"), true);
    program1ExtendSeeds.outAlignments = new PrintWriter(new FileWriter(arg[1] + ".align"), true);
    program1ExtendSeeds.outFragmentSizeAnalysis = new PrintWriter(new FileWriter(arg[1] + ".fragSizes"), true);
    program1ExtendSeeds.outOtherSolutions = new PrintWriter(new FileWriter(arg[1] + ".otherSolutions"), true);
   
   


    // Set constants:
    
    program1ExtendSeeds.program2DP.MAXTOTALJOINEDFRAGMENTSALLOWED = Fragment.MAXTOTALJOINEDFRAGMENTSALLOWED;
    program1ExtendSeeds.program2DP.maxNoOfConsecutiveFragmentsAlignedTogetherI = Fragment.maxNoOfConsecutiveFragmentsAlignedTogetherI;
    program1ExtendSeeds.program2DP.maxNoOfConsecutiveFragmentsAlignedTogetherJ = Fragment.maxNoOfConsecutiveFragmentsAlignedTogetherJ;


    program1ExtendSeeds.program2DP.SOMAV3 = false;    


    // Load reference maps:
    program1ExtendSeeds.program2DP.loadReferenceFragments(arg[2]);

    
    if (arg[3].equals("score"))
      ExtendSeeds.COMPARISON_BY_PVALUE = false;
    else if (arg[3].equals("pvalue"))
      ExtendSeeds.COMPARISON_BY_PVALUE = true;
    else
    { System.err.println("Please choose between 'pvalue' and 'score'!");
      System.exit(0);
    }
    
    // Load optical maps:
    if (arg[4].equals("allMaps"))
        program1ExtendSeeds.loadOMs(arg[0],-1,-1);
    else if (arg[4].equals("select"))
        program1ExtendSeeds.loadOMs(arg[0],Integer.parseInt(arg[5]),Integer.parseInt(arg[6]));
    else
    { System.err.println("Please choose between 'allMaps' and 'select'!");
      System.exit(0);
    }
    
    program1ExtendSeeds.TOTALFRAGMENTS = program1ExtendSeeds.program2DP.TOTALFRAGMENTS;

    // Aling the optical maps against the Reference genome:
    program1ExtendSeeds.alignOpticalMaps(0,program1ExtendSeeds.noOpticalMaps - 1);
    
    // Print statistics
    //program1ExtendSeeds.printStatistics();
    


    program1ExtendSeeds.out.println("\n\nTotal time:  " + (float)((double)(System.currentTimeMillis() - timeIn) / 3600000) + " h");
    
    program1ExtendSeeds.out.flush();
    program1ExtendSeeds.out.close();
    
    program1ExtendSeeds.outOK.flush();
    program1ExtendSeeds.outOK.close();
    
    program1ExtendSeeds.outNotOK.flush();
    program1ExtendSeeds.outNotOK.close();
    
    program1ExtendSeeds.outNotFound.flush();
    program1ExtendSeeds.outNotFound.close();
    
    program1ExtendSeeds.outAlignments.flush();
    program1ExtendSeeds.outAlignments.close();
    
    program1ExtendSeeds.outDiscarded.flush();
    program1ExtendSeeds.outDiscarded.close();
    
    program1ExtendSeeds.outFragmentSizeAnalysis.flush();
    program1ExtendSeeds.outFragmentSizeAnalysis.close();
    
    program1ExtendSeeds.outOtherSolutions.flush();
    program1ExtendSeeds.outOtherSolutions.close();
  }

}
