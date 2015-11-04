/*  Type: Fragment f-v. optimized for BioNano data

    Set main parameters and fragment information.

    
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

import java.util.ArrayList;
import java.util.List;
import java.util.Arrays;
import java.util.Comparator;
import java.text.DecimalFormat;


public class Fragment implements Comparable<Fragment>
{ 
  public static double CSIGMASEEDS = 1.2; //0.9 //2.0
  public static double CSIGMA = 2.6; //2.8 //3.0
  public static double BOUNDARYBIG = 1.9; //2.0
  public static int BOUNDARYBIGCONSTANT = 3;
  public static double BOUNDARYSMALL = 1.7; //2.0;
  public static int BOUNDARYSMALLCONSTANT = 3; //4;
  public static int MAXSIZESMALLFRAGMENTS = 2000; // //2050
  public static int MAXSIZEVERYSMALLFRAGMENTS = 800; // //880
  
  public static int LOWNUMBEROFFRAGMENTS = 10;
  public static int MINMOLECULESIZE = 50000; //150000

  public static int MAXALLOWEDMATCHESPERSEED = 70000; //30000
  public static int FALSECUTSSCORE = 1;
  public static int NUMBEROFFEASIBLESOLUTIONS = 6000000;  //30000
  public static int EXPECTEDNUMBEROF2SEEDMATCHES = 70000; //30000
  public static int MINIMUMNUMBEROFFEASIBLESOLUTIONS = 60;
  
  public static double ADJUSTMENT_OPTICAL_MAP_SIZE = 1.0;

  public static double MAXTOTALJOINEDFRAGMENTSALLOWED = 0.7;
  public static double MAXTOTALJOINEDFRAGMENTSALLOWED_STRINGENT = 0.7;
  public static int maxNoOfConsecutiveFragmentsAlignedTogetherI = 3; //5
  public static int maxNoOfConsecutiveFragmentsAlignedTogetherIFlexible = 3; //5
  public static int maxNoOfConsecutiveFragmentsAlignedTogetherJ = 5; //8
  public static double MAXTOTALEXTRACUTSALLOWED = 6;
  public static double MAXTOTALEXTRACUTSALLOWED_STRINGENT = 6;

  public static double CHISQUAREREWARDCONSTANT = 2;
  public static double CHISQUAREREWARDCONSTANTFORGAPS = 1.2;
  public static double SIMILARITYSCORES = 1 * Math.pow(CSIGMA, 2.0) * RestrictionMap.maxNoOfMatches;
  
  boolean fromOpticalMap = false;
  int size;
  int originalSize;
  int stDev = 0;
  int seedStDevLowerBound = 0;
  int seedStDevUpperBound = 0;
  int stDevLowerBound = 0;
  int stDevUpperBound = 0;
  int squaredStDev = 0;
  int location = -1; // in fragments
  long locationInBases = -1;
  boolean isInverted = false;  // Ref seed fragments
  ArrayList<Fragment> fragmentsMatched;
  ArrayList<Fragment> newFragmentsMatched;
  RestrictionMap opticalMap;
  
  public static double PVALUE = 2.24e-5;
  public static double DIFFERENCESECONDBESTPVALUE = 5.0;
  public static final int KMERSIZE = 3;  
  public static final int MAX_NUMBER_OTHER_SOLUTIONS = 10000000;
  public static final int FREQUENCY_OTHER_SOLUTIONS_TO_PLOT_FOR_QVALUE_ANALYSIS = 20;
  public static final DecimalFormat NUMBER_FORMAT = new DecimalFormat("#.000");  


  public Fragment()  // Empty fragment
  {  size = 0;
  }

  public Fragment(int aSize, RestrictionMap anOpticalMap, int aLocation)  // OM
  {  size = aSize;
     if (size < 0) size = 0;
     stDev = 0;
     location = aLocation;
     fromOpticalMap = true;
     opticalMap = anOpticalMap;
  }
  
  public Fragment(int aSize, int aLocation, long aLocationInBases, boolean fragIsInverted)   // Reference maps fragments
  {  location = aLocation;
     size = aSize;
     stDev = chooseBinStDev(aSize);
     locationInBases = aLocationInBases;
     isInverted = fragIsInverted;
     originalSize = aSize;

     // For OpGen's data:
     /*
     if (aSize <= 20000)
     { size = (int)(aSize * 0.975 + 100);
     }
     else
     { size = (int)(aSize - 400);     
     }
     */
     
     
     squaredStDev = (int)Math.pow((double)stDev, 2.0);
     stDevUpperBound = size - (int)(CSIGMA * stDev);
     stDevLowerBound = size + (int)(CSIGMA * stDev);
     seedStDevUpperBound = size - (int)(CSIGMASEEDS * stDev);
     seedStDevLowerBound = size + (int)(CSIGMASEEDS * stDev);
  }
 
  public int compareTo(Fragment y) // y is a fragment from reference genome (i.e., with stdev != 0)
  { if (this.size < y.size)
    { return -1; }
    else if (this.size > y.size)
    { return 1; }
    return 0;
  }

  public int compareToCSIGMA(Fragment y) // y is a fragment from reference genome (i.e., with stdev != 0)
  { if (this.size < y.stDevLowerBound)
    { return -1; }
    else if (this.size > y.stDevUpperBound)
    { return 1; }
    return 0;
  }


  public int chooseBinStDev(int aSize)
  { // For OpGen's data:  
    // return (int)(0.03*aSize + 450);
    
    // For BioNano Genomics' data:
       return (int)(0.02*aSize + 200);

    /* // Default:
    if (aSize < 800)
    { return (int) (aSize*0.25 + 200);
    }
    if (aSize < 2000)
    { return (int) (aSize*0.41);
    }
    if (aSize < 4000)
    { return (int) (aSize*0.11);
    }
    if (aSize < 8000)
    { return (int) (aSize*0.095);
    }
    if (aSize < 20000)
    { return (int) (aSize*0.08);
    }
    // aSize > 20000
    return (int) (aSize*0.0625);
    */
  }


  public int strictCompareToSeedsUpperBound(Fragment y) //, Fragment h) // y is a Fragment from reference genome (i.e., with stdev != 0)  // h next one
  { if (this.size < y.seedStDevUpperBound)
    { return -1; }
    if (this.size > y.seedStDevUpperBound)
    { return 1; }

    return 0;
  }

  public int strictCompareToSeedsLowerBound(Fragment y) //, Fragment h) // y is a Fragment from reference genome (i.e., with stdev != 0)   // h previous one
  { if (this.size < y.seedStDevLowerBound)
    { return -1; }
    if (this.size > y.seedStDevLowerBound)
    { return 1; }

    return 0;
  }
  
  public int strictCompareToUpperBound(Fragment y, int deviation) // y is a Fragment from reference genome (i.e., with stdev != 0)  // h next one  Fragment h,
  { if (this.size < (y.size - deviation))
    { return -1; }
    if (this.size > (y.size - deviation))
    { return 1; }

    return 0;
  }
  
  public int strictCompareToLowerBound(Fragment y, int deviation) // y is a Fragment from reference genome (i.e., with stdev != 0)   // h previous one  Fragment h,
  { if (this.size < (y.size + deviation))
    { return -1; }
    if (this.size > (y.size + deviation))
    { return 1; }

    return 0;
  }


  public String toString()
  { return "" + originalSize;
  }
   


  // Advanced methods:

  public ArrayList<Fragment> binarySearchWithList(ArrayList<Fragment> listReferenceSeeds, Fragment[] listReferenceSeedsToArray, int low, int high) // y is a Fragment from OM
  { // return the number of matches of y in the list   
    int highest = high;  // exact index of latest value
    ArrayList<Fragment> fragmentsMatched = new ArrayList<Fragment>();

    int lowerBoundList = justBinarySearch(listReferenceSeeds, low, high, false);
    if (lowerBoundList == -1)
    { return fragmentsMatched;
    }

    int upperBoundList = justBinarySearch(listReferenceSeeds, low, high, true);
    if (upperBoundList == -1)
    { return fragmentsMatched;
    }
    
    if (upperBoundList - lowerBoundList < 0)
    { return fragmentsMatched;
    }

    Fragment[] tmpFragments2 = new Fragment[upperBoundList - lowerBoundList + 1];
    System.arraycopy(listReferenceSeedsToArray, lowerBoundList, tmpFragments2, 0, upperBoundList - lowerBoundList + 1);
    List <Fragment> list = Arrays.asList(tmpFragments2);
    fragmentsMatched = new ArrayList<Fragment>(upperBoundList - lowerBoundList + 1);
    fragmentsMatched.addAll(list);

    return fragmentsMatched;
  }



  public int justBinarySearch(ArrayList<Fragment> listReferenceSeeds, int low, int high, boolean upperBound) // y is a Fragment from OM
  { // return a location in listReferenceSeeds given y
    Fragment y = this;

    int highest = high;  // exact index of latest value
    int mid = -1;

    // continue searching while [imin,imax] is not empty
    while (high >= low)
    { mid = low + (high - low) / 2;
      int result = 0; 

      if (! upperBound)
      { result = y.strictCompareToSeedsLowerBound(listReferenceSeeds.get(mid));
      }
      else
      { result = y.strictCompareToSeedsUpperBound(listReferenceSeeds.get(mid)); //, h);
      }

      if (result < 0)
      { high = mid - 1; }
      else if (result > 0)
      { low = mid + 1; }
      else
      { low = mid;
        high = mid - 1;
      }
    }
    
    if ( upperBound )
    { if (y.strictCompareToSeedsUpperBound(listReferenceSeeds.get(mid)) < 0)
      { if (mid > 0) mid--; }
      
      while (mid < highest && y.strictCompareToSeedsUpperBound(listReferenceSeeds.get(mid+1)) >= 0)
      { mid++;
      }
      
      if (y.strictCompareToSeedsUpperBound(listReferenceSeeds.get(mid)) < 0)
      { return -1;
      }
    }
    else
    { if (y.strictCompareToSeedsLowerBound(listReferenceSeeds.get(mid)) > 0)
      { if (mid < highest) mid++; }
      
      while (mid > 0 && y.strictCompareToSeedsLowerBound(listReferenceSeeds.get(mid-1)) <= 0)
      { mid--;
      }
      
      if (y.strictCompareToSeedsLowerBound(listReferenceSeeds.get(mid)) > 0)
      { return -1;
      }
    }

    return mid;
  }
  

  
  // override
  public int justBinarySearch(ArrayList<Fragment> listReferenceSeeds, int low, int high, boolean upperBound, int deviation) // y is a Kmer from OM
  { // return a location in listReferenceSeeds given y
    Fragment y = this;

    int highest = high;  // exact index of latest value
    int mid = -1;

    // continue searching while [imin,imax] is not empty
    while (high >= low)
    { mid = low + (high - low) / 2;
      int result = 0;
      if (! upperBound)
      { result = y.strictCompareToLowerBound(listReferenceSeeds.get(mid), deviation); //h,
      }
      else
      { result = y.strictCompareToUpperBound(listReferenceSeeds.get(mid), deviation);  // h,
      }

      if (result < 0)
      { high = mid - 1; }
      else if (result > 0)
      { low = mid + 1; }
      else
      { low = mid;
        high = mid - 1;
      }
    }
    
    if ( upperBound )
    { if (y.strictCompareToUpperBound(listReferenceSeeds.get(mid), deviation) < 0)
      { if (mid > 0) mid--; }
      
      while (mid < highest && y.strictCompareToUpperBound(listReferenceSeeds.get(mid+1), deviation) >= 0)
      { mid++;
      }
      
      if (y.strictCompareToUpperBound(listReferenceSeeds.get(mid), deviation) < 0)
      { return -1; }
    }
    else
    { if (y.strictCompareToLowerBound(listReferenceSeeds.get(mid), deviation) > 0)
      { if (mid < highest) mid++; }
      
      while (mid > 0 && y.strictCompareToLowerBound(listReferenceSeeds.get(mid-1), deviation) <= 0)
      { mid--;
      }
      
      if (y.strictCompareToLowerBound(listReferenceSeeds.get(mid), deviation) > 0)
      { return -1; }
    }

    return mid;
  }
  


  public static class Comparators
  {
    public static Comparator<Fragment> MATCHINGLISTSIZE = new Comparator<Fragment>() {
        @Override
        public int compare(Fragment frag1, Fragment frag2) {
            if (frag1.newFragmentsMatched.size() < frag2.newFragmentsMatched.size())
              return -1;
            if (frag1.newFragmentsMatched.size() > frag2.newFragmentsMatched.size())
              return 1;
            return 0;
        }
    };
  }

}
