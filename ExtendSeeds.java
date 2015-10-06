/*  ExtendSeeds v8

    Align each optical map to the reference maps by using a seeding approach
    together with a technology-independent statistical analysis.
    
    
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

import java.util.Scanner;
import java.io.File;
import java.io.IOException;
import java.io.Writer;
import java.io.PrintWriter;
import java.io.FileWriter;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Arrays;
import java.util.Random;

import cern.jet.random.Binomial;
import cern.jet.random.Poisson;
import cern.jet.random.Normal;
import cern.jet.random.ChiSquare;
import cern.jet.random.engine.RandomEngine;
import cern.jet.random.engine.DRand;

import cern.jet.stat.Descriptive;
import cern.jet.random.Normal;
import cern.colt.list.DoubleArrayList;
import cern.colt.list.IntArrayList;


public class ExtendSeeds
{
  public class CandidateSolution implements Comparable<CandidateSolution>
  { 
    int index;

    double score = RestrictionMap.MINIMUMSCORE;

    int location = -1;
    boolean isInvertedStrand = false;

    Fragment yOM = null;
    Fragment zRef = null;

    int lengthInReference = 0;
    int mapMatches = 0;
    int refSizeMinusMissingFragments = 0;
    int noSmallRefFragments = 0;
    
    boolean isUnique = true;

    int missingCuts = 0;
    int falseCuts = 0;
    int missingFragments = 0;
    double chiSquare = 0.0;
    double chiSquarePerMatch = 0.0;   // normalized by Wilson-Hilferty Transformation

    int usefulRefSize = 0;
    int usefulOMSize = 0;
    double usefulSizeComparison = 0.0;

    double probMissingCuts = 0.0;
    double probFalseCuts = 0.0;
    double probChiSquarePerMatch = 0.0;
      
    double zScore = 0.0;    // Total Z-Score
    double zScoreEC = 0.0;  // Error cuts (missing cuts + false cuts + missing fragments)
    double zScoreFCMM = 0.0; // False cuts / Inverse(Map Matches)
    double zScoreCSPM = 0.0; // Chi Square per match (normalized by Wilson-Hilferty Transformation)
    double pValue = 0.0;  // Computed from Z-Score
    
    double secondBestPValue = 0.0;
    int secondBestLocation = -1;

    double avgFragmentSizeOM = 0.0;
    double avgFragmentSizeOMMinusEnds = 0.0;
    double avgFragmentSizeReferenceMinusEnds = 0.0;
    double fragSizeErrorAbs = 0.0;
    double fragSizeErrorRelative = 0.0;

    int countFragmentsOM = 0;
    int countFragmentsOMMinusEnds = 0;
    int countFragmentsReference = 0;
    int countFragmentsReferenceMinusEnds = 0;
    
    int verySmallFragmentsNotDetectedAsMissing = 0;
    int verySmallFragmentsDetectedAsMissing = 0;

    int seedsInSameLocation = 0;

    boolean significant = false;
    
    RestrictionMap opticalMap = null;
    
    public CandidateSolution(int index, RestrictionMap opticalMap)
    {  this.index = index;
       this.opticalMap = opticalMap;
    }


    public CandidateSolution(int index, RestrictionMap opticalMap, double pValue)
    {  this.index = index;
       this.opticalMap = opticalMap;
       this.pValue = pValue;
    }


    public int compareTo(CandidateSolution obj2)
    {  
       if (COMPARISON_BY_PVALUE)
       { if (this.pValue < obj2.pValue)
           return -1;
         if (this.pValue > obj2.pValue)
           return 1;
         return 0;
       }
       else
       { if (this.score > obj2.score)
           return -1;
         if (this.score < obj2.score)
           return 1;
         return 0; 
       }
    }

  }



  boolean PRINTFRAGMENTS = false;

  int TOTALFRAGMENTS;
  public static boolean COMPARISON_BY_PVALUE; //final

  int alignedMaps = 0;
  int noInvertedMaps = 0;
  int discardedMaps = 0;
  int noAlignmentFound = 0;
  int nonUniqueAlignmentFound = 0;
  
  CandidateSolution bestAlignment = null;



  PrintWriter out, outOK, outNotOK, outNotFound, outDiscarded, outAlignments, outOtherSolutions, outFragmentSizeAnalysis;

  cern.jet.random.engine.DRand generator = new cern.jet.random.engine.DRand();
  cern.jet.random.ChiSquare chiSquareDistribution = new cern.jet.random.ChiSquare(6,generator);


  RestrictionMap[] opticalMaps;
  int noOpticalMaps = 0;
  int originalNoOpticalMaps = 0;
  RestrictionMap opticalMap;
  
  int totalNoOfMaps = 0;
  int noNotAlignedMaps = 0;
  int noAlignedMaps = 0;
  int notUniqueAlignedMaps = 0;
  int noAlignedMapsWithoutSignificance = 0;
  
  long totalSizeMaps = 0;

  int matchesSelectedFragment = 0;
  boolean fragmentMatchesAlreadyPrinted = false;
  int noFeasibleAlignmentsCurrentOMFragment = 0;

  int feasibleAlignments = 0;
  int seedsTried = 0;

  String headerAlignment = "";

  cern.jet.random.Binomial binomialDistributionMissingCuts;
  
  cern.jet.random.Poisson poissonDistributionExtraCuts;
  double mean = 0; //Fragment.ROOTCHISQUAREDISTRIBUTIONMEAN; 
  double sd = 1; //Fragment.ROOTCHISQUAREDISTRIBUTIONSD;
  
  cern.jet.random.Normal normalDistributionChiSquarePerMatch = new cern.jet.random.Normal(mean, sd, generator);
  cern.jet.random.Normal standardNormalDistribution = new cern.jet.random.Normal(0, 1, generator);

  org.apache.commons.math3.distribution.ExponentialDistribution exponentialDistributionChiSquarePerMatch = new org.apache.commons.math3.distribution.ExponentialDistribution(2);


  DynamicProgramming program2DP = new DynamicProgramming();
  
  
  volatile int countOtherSolutionsPlot = 0;
  
  
  // For SOMA-like score:
  boolean errorOnMAXIT = false;
  static final double FDISTRIBPVALUE = 0.05; //0.01



  public static void main (String[] arg) throws IOException
  { if (arg.length < 6)
    { System.err.println("Use: java FullAlignment/ExtendSeeds Optical.maps(Optical Maps) Output.txt(fileNameOutput) locationInSilico_IfHuman/InSilicoFile_Otherwise human? allMaps? [pvalue|score] (OPTIONAL firstIndex lastIndex)");
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
   
   


    // CONSTANTS:
    
    program1ExtendSeeds.program2DP.MAXTOTALJOINEDFRAGMENTSALLOWED = Fragment.MAXTOTALJOINEDFRAGMENTSALLOWED;
    program1ExtendSeeds.program2DP.maxNoOfConsecutiveFragmentsAlignedTogetherI = Fragment.maxNoOfConsecutiveFragmentsAlignedTogetherI;
    program1ExtendSeeds.program2DP.maxNoOfConsecutiveFragmentsAlignedTogetherJ = Fragment.maxNoOfConsecutiveFragmentsAlignedTogetherJ;


    program1ExtendSeeds.program2DP.SOMAV3 = false;    

    // Load Reference genome:
    if ( arg[3].equals("true") )
      program1ExtendSeeds.program2DP.loadReferenceFragmentsHuman(arg[2]);
    else
      program1ExtendSeeds.program2DP.loadReferenceFragments(arg[2]);

    // Load optical maps:    
    if (arg[5].equals("score"))
      COMPARISON_BY_PVALUE = false;
    else if (arg[5].equals("pvalue"))
      COMPARISON_BY_PVALUE = true;
    else
    { System.err.println("Please choose between 'pvalue' and 'score'!");
      System.err.println("Use: java FullAlignment/ExtendSeeds Optical.maps(Optical Maps) Output.txt(fileNameOutput) locationInSilico_IfHuman/InSilicoFile_Otherwise human? allMaps? [pvalue|score] (OPTIONAL firstIndex lastIndex)");
      System.exit(0);
    }

    
    if ( Boolean.parseBoolean(arg[4]) )
        program1ExtendSeeds.loadOMs(arg[0],-1,-1);
    else
        program1ExtendSeeds.loadOMs(arg[0],Integer.parseInt(arg[6]),Integer.parseInt(arg[7]));

    program1ExtendSeeds.TOTALFRAGMENTS = program1ExtendSeeds.program2DP.TOTALFRAGMENTS;

    // Aling the optical maps against the Reference genome:
    program1ExtendSeeds.alignOpticalMaps(0,program1ExtendSeeds.noOpticalMaps - 1);


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



  public void alignOpticalMaps(int FIRSTOMTOALIGN, int LASTOMTOALIGN)
  { for (int noOM = FIRSTOMTOALIGN; noOM <= LASTOMTOALIGN; noOM++)
    { RestrictionMap opticalMap = opticalMaps[noOM];

      PRINTFRAGMENTS = false;
      totalNoOfMaps++;

      boolean aligned = alignMap(opticalMap);
    }

  }



  public boolean alignMap(RestrictionMap opticalMap)
  {
    bestAlignment = null;
    
    feasibleAlignments = 0;
    seedsTried = 0;
    
    ArrayList<CandidateSolution> listFeasibleAlignments = new ArrayList<CandidateSolution>(Fragment.NUMBEROFFEASIBLESOLUTIONS);

    boolean firstPrintPerMap = true;  // print heading for the map and its number of seed matches

    headerAlignment = "";
    int i = 0;
    int j = 0;

    // Find seeds:
    
    int highestOMFragmentToAnalyze = (opticalMap.size() - 2);   // do not consider OM ending fragments

    int shortCount = 0;
    ArrayList<Fragment> listAllFragments = new ArrayList<Fragment>(opticalMap.noFragments - 2);

    for (i = 0; i < opticalMap.size(); i++)   
    { Fragment y = opticalMap.get(i);
    
      // do not consider OM ending fragments
      if (i == 0 || i > highestOMFragmentToAnalyze || shortCount > 15) // speed-up
      { y.newFragmentsMatched = new ArrayList<Fragment>(1);
        //shortCount++; // !!!
        continue;
      }
      shortCount++;
      
      
      y.newFragmentsMatched = new ArrayList<Fragment>(Fragment.EXPECTEDNUMBEROF2SEEDMATCHES);
      ArrayList<Fragment> fragmentsMatched = y.binarySearchWithList(program2DP.listReferenceSeeds, program2DP.listReferenceSeedsToArray, 0, (program2DP.noReferenceSeeds - 1));
      int matches = fragmentsMatched.size();
      
      if (matches > 0) // && matches <= Fragment.MAXALLOWEDMATCHESPERSEED)
      { 
          y.newFragmentsMatched = fragmentsMatched;
          listAllFragments.add(y);
      }
    } 
    
    
YFor: for (Fragment y : listAllFragments)
    {
      ArrayList<Fragment> fragmentsMatched = y.newFragmentsMatched;
      i = y.location;

      int matches = fragmentsMatched.size();

      if (matches == 0)
      { continue YFor;
      }

      int dimensionsLeft = i;
      int dimensionsRight = opticalMap.size() - (i + 1);
      Fragment[] opticalMapLeft = new Fragment[dimensionsLeft];
      Fragment[] opticalMapRight = new Fragment[dimensionsRight];

      // Both Forward and Reverse strands, in case the enzyme is palindrome:
      System.arraycopy(opticalMap.fragments, 0, opticalMapLeft, 0, dimensionsLeft);
      List <Fragment> list = Arrays.asList(opticalMapLeft);
      Collections.reverse(list);
      opticalMapLeft = (Fragment[]) list.toArray();
      System.arraycopy(opticalMap.fragments, i + 1, opticalMapRight, 0, dimensionsRight);  // + 1

      fragmentMatchesAlreadyPrinted = false;    // for list of seed matches

      noFeasibleAlignmentsCurrentOMFragment = 0;

           

ZFor: for (Fragment z : fragmentsMatched) // Seeds matched in reference genome
      { 
        // Check if the OM is near the borders of the genome
        if (! z.isInverted)
        { if (((z.location == TOTALFRAGMENTS - 1) && (dimensionsRight > 0)) || ((z.location + dimensionsRight) > TOTALFRAGMENTS - 1 + Fragment.maxNoOfConsecutiveFragmentsAlignedTogetherIFlexible + 3))
          { continue ZFor; }
          if (((z.location == 0) && (dimensionsLeft > 0)) || (z.location - dimensionsLeft < 0 - (Fragment.maxNoOfConsecutiveFragmentsAlignedTogetherIFlexible + 3)))
          { continue ZFor; }
        }
        else
        { if (((z.location == TOTALFRAGMENTS - 1) && (dimensionsLeft > 0)) || ((z.location + dimensionsLeft) > TOTALFRAGMENTS - 1 + Fragment.maxNoOfConsecutiveFragmentsAlignedTogetherIFlexible + 3))
          { continue ZFor; }
          if (((z.location == 0) && (dimensionsRight > 0)) || (z.location - dimensionsRight < 0 - (Fragment.maxNoOfConsecutiveFragmentsAlignedTogetherIFlexible + 3)))
          { continue ZFor; }
        }
        
        seedsTried++;

        // Compute score for local seed match:
        double chiSquareSeedMatchAbs = Math.abs(y.size - z.size) / (double)z.stDev;
        double chiSquareSeedMatch = Math.pow(chiSquareSeedMatchAbs, 2.0);

        double newScore = -chiSquareSeedMatch;
        DynamicProgramming.Pair pLeft = program2DP.new Pair(0, -1);
        DynamicProgramming.Pair pRight = program2DP.new Pair(0, -1);

        int candidateMissingCuts = 0;
        int candidateFalseCuts = 0;
        int candidateMissingFragments = 0;

        int currentUsefulRefSize = 0;
        int currentUsefulOMSize = 0;

        int currentRefSizeMinusMissingFragments = 0;
        
        int currentNoSmallRefFragments = 0;

        double averagedMissingCutsI = 0;
        double averagedMissingCutsJ = 0;
        int lengthJ = 0;
        
        
        double avgFragmentSizeOM = 0.0;
        double avgFragmentSizeOMMinusEnds = 0.0;
        double avgFragmentSizeReferenceMinusEnds = 0.0;
        double fragSizeErrorAbs = 0.0;
        double fragSizeErrorRelative = 0.0;
        
        int countFragmentsOM = 0;
        int countFragmentsOMMinusEnds = 0;
        int countFragmentsReference = 0;
        int countFragmentsReferenceMinusEnds = 0;
        
        int verySmallFragmentsNotDetectedAsMissing = 0;
        int verySmallFragmentsDetectedAsMissing = 0;
        
        if (dimensionsLeft > 0)
        { pLeft = program2DP.dynamicProgrammingSOMAv2(opticalMapLeft, ( z.isInverted ? z.location + 1 : z.location - 1), (! z.isInverted), TOTALFRAGMENTS - 1);  // third parameter should be true in case of a non-inverted match and false in case of an inverted match
          if (pLeft.bestScore <= RestrictionMap.MINIMUMSCORE)
          { continue ZFor;
          }

          newScore += pLeft.bestScore;

          candidateMissingCuts = program2DP.currentMissingCuts;
          candidateFalseCuts = program2DP.currentFalseCuts;
          candidateMissingFragments = program2DP.currentMissingFragments;
          
          currentNoSmallRefFragments = program2DP.currentNoSmallRefFragments;

          lengthJ = program2DP.locationMaxScoreJ;

          currentRefSizeMinusMissingFragments = program2DP.currentRefSizeMinusMissingFragments;

          currentUsefulRefSize = program2DP.currentUsefulRefSize;
          currentUsefulOMSize = program2DP.currentUsefulOMSize;
          
          avgFragmentSizeOM += program2DP.avgFragmentSizeOM;
          avgFragmentSizeOMMinusEnds += program2DP.avgFragmentSizeOMMinusEnds;
          avgFragmentSizeReferenceMinusEnds += program2DP.avgFragmentSizeReferenceMinusEnds;
          fragSizeErrorAbs += program2DP.fragSizeErrorAbs;
          fragSizeErrorRelative += program2DP.fragSizeErrorRelative;
        
          countFragmentsOM += program2DP.countFragmentsOM;
          countFragmentsOMMinusEnds += program2DP.countFragmentsOMMinusEnds;
          countFragmentsReference += program2DP.countFragmentsReference;
          countFragmentsReferenceMinusEnds += program2DP.countFragmentsReferenceMinusEnds;
          
          verySmallFragmentsNotDetectedAsMissing += program2DP.verySmallFragmentsNotDetectedAsMissing;
          verySmallFragmentsDetectedAsMissing += program2DP.verySmallFragmentsDetectedAsMissing;
        }

        if (dimensionsRight > 0)
        { pRight = program2DP.dynamicProgrammingSOMAv2(opticalMapRight, ( z.isInverted ? z.location - 1 : z.location + 1), z.isInverted, TOTALFRAGMENTS - 1);  // third parameter should be false in case of a non-inverted match and true in case of an inverted match
          if (pRight.bestScore <= RestrictionMap.MINIMUMSCORE)
          { continue ZFor; }

          newScore += pRight.bestScore;
          
          candidateMissingCuts += program2DP.currentMissingCuts;
          candidateFalseCuts += program2DP.currentFalseCuts;
          candidateMissingFragments += program2DP.currentMissingFragments;
          
          currentNoSmallRefFragments += program2DP.currentNoSmallRefFragments;
          
          lengthJ += program2DP.locationMaxScoreJ;

          currentRefSizeMinusMissingFragments += program2DP.currentRefSizeMinusMissingFragments;

          currentUsefulRefSize += program2DP.currentUsefulRefSize;
          currentUsefulOMSize += program2DP.currentUsefulOMSize;
          
          avgFragmentSizeOM += program2DP.avgFragmentSizeOM;
          avgFragmentSizeOMMinusEnds += program2DP.avgFragmentSizeOMMinusEnds;
          avgFragmentSizeReferenceMinusEnds += program2DP.avgFragmentSizeReferenceMinusEnds;
          fragSizeErrorAbs += program2DP.fragSizeErrorAbs;
          fragSizeErrorRelative += program2DP.fragSizeErrorRelative;
        
          countFragmentsOM += program2DP.countFragmentsOM;
          countFragmentsOMMinusEnds += program2DP.countFragmentsOMMinusEnds;
          countFragmentsReference += program2DP.countFragmentsReference;
          countFragmentsReferenceMinusEnds += program2DP.countFragmentsReferenceMinusEnds;
          
          verySmallFragmentsNotDetectedAsMissing += program2DP.verySmallFragmentsNotDetectedAsMissing;
          verySmallFragmentsDetectedAsMissing += program2DP.verySmallFragmentsDetectedAsMissing;
        }

        lengthJ++;
        
        currentRefSizeMinusMissingFragments += z.size;
        
        currentUsefulRefSize += z.size;
        currentUsefulOMSize += y.size;
        
        avgFragmentSizeOM += y.size;
        avgFragmentSizeOMMinusEnds += y.size;
        avgFragmentSizeReferenceMinusEnds += z.size;
        fragSizeErrorAbs += (Math.abs(y.size - z.size) / (double)z.size);
        fragSizeErrorRelative += ((y.size - z.size) / (double)z.size);

        countFragmentsOM++;
        countFragmentsOMMinusEnds++;
        countFragmentsReference++;
        countFragmentsReferenceMinusEnds++;




        if (newScore <= RestrictionMap.MINIMUMSCORE)
        { continue ZFor;
        }
        

        // We have found a candidate solution!  Now compute its 'forward' location and number of fragments covered in the reference maps:
        int candidateLocation = -1;
        int candidateLengthInReference = 0;
        if (! z.isInverted)
        { candidateLocation = (dimensionsLeft > 0 ? pLeft.bestLastLocationJ : z.location);
          candidateLengthInReference = (dimensionsRight > 0 ? pRight.bestLastLocationJ : z.location) - candidateLocation + 1;
        }
        else
        { candidateLocation = (dimensionsRight > 0 ? pRight.bestLastLocationJ : z.location);
          candidateLengthInReference = (dimensionsLeft > 0 ? pLeft.bestLastLocationJ : z.location) - candidateLocation + 1;
        }
        

        // Statistical significance tests:
        int currentMapMatches = opticalMap.size() - candidateFalseCuts;
        
        double currentChiSquare = (Math.abs(newScore + program2DP.matchBonus * (Fragment.FALSECUTSSCORE*candidateFalseCuts + candidateMissingCuts + candidateMissingFragments + verySmallFragmentsNotDetectedAsMissing + verySmallFragmentsDetectedAsMissing)));
        
        double currentChiSquarePerMatchWilsonHilfertyTransformation = (Math.cbrt(currentChiSquare / (currentMapMatches-2)) - (1 - 2.0/(9*(currentMapMatches-2)))) / Math.sqrt(2.0 / (9 * (currentMapMatches-2)));
        double currentChiSquarePerMatch = currentChiSquarePerMatchWilsonHilfertyTransformation;
        
        
              
        // First checks:
        if ((candidateMissingCuts / (double)(candidateLengthInReference - candidateMissingFragments - verySmallFragmentsNotDetectedAsMissing - verySmallFragmentsDetectedAsMissing - 1.0)) > (Fragment.MAXTOTALJOINEDFRAGMENTSALLOWED) || (candidateFalseCuts / (currentRefSizeMinusMissingFragments/100000.0) > Fragment.MAXTOTALEXTRACUTSALLOWED))
          continue ZFor;        
        
        
        
        // Store solution
        
        //---------------------------

        
        CandidateSolution properties = new CandidateSolution(feasibleAlignments, opticalMap);
        feasibleAlignments++;
        listFeasibleAlignments.add(properties);

        properties.score = newScore;
        
        properties.lengthInReference = candidateLengthInReference; 
        properties.location = candidateLocation;
        properties.isInvertedStrand = z.isInverted;

        properties.missingCuts = candidateMissingCuts;
        properties.falseCuts = candidateFalseCuts;
        properties.missingFragments = candidateMissingFragments;
        
        properties.chiSquarePerMatch = currentChiSquarePerMatch;
        properties.noSmallRefFragments = currentNoSmallRefFragments;

        properties.mapMatches = currentMapMatches;

        properties.refSizeMinusMissingFragments = currentRefSizeMinusMissingFragments;

        properties.usefulRefSize = currentUsefulRefSize;
        properties.usefulOMSize = currentUsefulOMSize;
        
        properties.chiSquare = currentChiSquare;
        properties.yOM = y;
        properties.zRef = z;
        
        
        properties.avgFragmentSizeOM = avgFragmentSizeOM;
        properties.avgFragmentSizeOMMinusEnds = avgFragmentSizeOMMinusEnds;
        properties.avgFragmentSizeReferenceMinusEnds = avgFragmentSizeReferenceMinusEnds;
        properties.fragSizeErrorAbs = fragSizeErrorAbs;
        properties.fragSizeErrorRelative = fragSizeErrorRelative;
        
        properties.countFragmentsOM = countFragmentsOM;
        properties.countFragmentsOMMinusEnds = countFragmentsOMMinusEnds;
        properties.countFragmentsReference = countFragmentsReference;
        properties.countFragmentsReferenceMinusEnds = countFragmentsReferenceMinusEnds;
        
        properties.verySmallFragmentsNotDetectedAsMissing = verySmallFragmentsNotDetectedAsMissing;
        properties.verySmallFragmentsDetectedAsMissing = verySmallFragmentsDetectedAsMissing;
        
        //---------------------------
        
        
        noFeasibleAlignmentsCurrentOMFragment++;

      } // END FOR - Ref Seed Z

    } // END FOR - OM fragment Y

    
    


    // If no alignment has been found:
    if (feasibleAlignments == 0)
    { noAlignmentFound++;
      noNotAlignedMaps++;
      
      plotNotFoundSolution(opticalMap,false,feasibleAlignments,seedsTried);

      return false;
    }


    
    DoubleArrayList m1 = new DoubleArrayList(listFeasibleAlignments.size());
    DoubleArrayList m2 = new DoubleArrayList(listFeasibleAlignments.size());
    DoubleArrayList m3 = new DoubleArrayList(listFeasibleAlignments.size());
    DoubleArrayList q = new DoubleArrayList(listFeasibleAlignments.size());

    for (int w = 0; w < (listFeasibleAlignments.size()); w++)
    {  m2.add (listFeasibleAlignments.get(w).falseCuts);
       m3.add (listFeasibleAlignments.get(w).missingCuts + listFeasibleAlignments.get(w).missingFragments + listFeasibleAlignments.get(w).falseCuts);
       q.add(listFeasibleAlignments.get(w).chiSquarePerMatch); // normalized W-H!!    Math.abs( + 0.5)  i.e. - -0.50

    } 

    double m2Mean = cern.jet.stat.Descriptive.mean(m2);
    cern.jet.stat.Descriptive.standardize(m2, m2Mean, cern.jet.stat.Descriptive.standardDeviation(cern.jet.stat.Descriptive.sampleVariance(m2, m2Mean)));
    double m3Mean = cern.jet.stat.Descriptive.mean(m3);
    cern.jet.stat.Descriptive.standardize(m3, m3Mean, cern.jet.stat.Descriptive.standardDeviation(cern.jet.stat.Descriptive.sampleVariance(m3, m3Mean)));
    double qMean = cern.jet.stat.Descriptive.mean(q);
    cern.jet.stat.Descriptive.standardize(q, qMean, cern.jet.stat.Descriptive.standardDeviation(cern.jet.stat.Descriptive.sampleVariance(q, qMean)));
  


    double persM1 = 0; //pers.correlation(m1.elements(), q.elements()); // input: 2 double[]

    for (int w = 0; w < (listFeasibleAlignments.size()); w++)
    { m1.add(((m2.get(w)) + m3.get(w) + q.get(w)));
    }
    
    double m1Mean = cern.jet.stat.Descriptive.mean(m1);
    cern.jet.stat.Descriptive.standardize(m1, m1Mean, cern.jet.stat.Descriptive.standardDeviation(cern.jet.stat.Descriptive.sampleVariance(m1, m1Mean)));
    
    
   
   
       
    
    double bestPValue = 1;
    double bestZScore = -4.0;
    int bestIndex = -1;
    double secondBestPValue = 1;
    double secondBestZScore = -4.0;
    int secondBestIndex = -1;
    
    // Plot one solution with a certain frequency
	int currentCount = 1;

    
    for (int w = 0; w < (listFeasibleAlignments.size()); w++)
    { 
      double zScore = m1.get(w);
      CandidateSolution properties = listFeasibleAlignments.get(w);
      
      // First checks:
      
      if (Double.isNaN(zScore))
	    continue;
	    
	    
	  if (currentCount >= Fragment.FREQUENCY_OTHER_SOLUTIONS_TO_PLOT_FOR_QVALUE_ANALYSIS && this.countOtherSolutionsPlot <= Fragment.MAX_NUMBER_OTHER_SOLUTIONS){
		countOtherSolutionsPlot++; 
		outOtherSolutions.println(Fragment.NUMBER_FORMAT.format(zScore) + "\t" + opticalMap.size());
		currentCount = 0;
	  }
	  currentCount++;
			
      
      
      if ((COMPARISON_BY_PVALUE ? (zScore > secondBestZScore) : ((secondBestIndex != -1 && properties.compareTo(listFeasibleAlignments.get(secondBestIndex)) > 0))))  continue;
        
      double candidateScore = properties.score;

      int candidateLengthInReference = properties.lengthInReference;
      int candidateLocation = properties.location;
      boolean candidateIsInvertedStrand = properties.isInvertedStrand;
      int candidateMissingCuts = properties.missingCuts;
      int candidateFalseCuts = properties.falseCuts;
      int candidateMissingFragments = properties.missingFragments;
        
      double candidateChiSquarePerMatch = properties.chiSquarePerMatch;

      double candidateMapMatches = properties.mapMatches;

      int candidateRefSizeMinusMissingFragments = properties.refSizeMinusMissingFragments;
      
      int candidateVerySmallFragmentsNotDetectedAsMissing = properties.verySmallFragmentsNotDetectedAsMissing;
      int candidateVerySmallFragmentsDetectedAsMissing = properties.verySmallFragmentsDetectedAsMissing;

      
      // M: Ex first checks:
               
      
      // Check map size:      
      int usefulRefSize = properties.usefulRefSize;
      int usefulOMSize = properties.usefulOMSize;


      double usefulSizeComparison = usefulOMSize / (double)usefulRefSize;
      
      // M: Useful size comparison
      
      properties.usefulSizeComparison = usefulSizeComparison;
  
      // M: Statistical significance tests     
           
	  double probabilityErrorsMC = 0;	
	  double probabilityErrorsFC = 0;
	  double pdfCSPMWH = 0;
	  double probabilityErrorsES = 0;
	  
	  double currentPValue = standardNormalDistribution.cdf(zScore);  

     
      properties.zScore = zScore;
      properties.zScoreEC = m3.get(w);
      properties.zScoreFCMM = m2.get(w);
      properties.zScoreCSPM = q.get(w);
      properties.pValue = currentPValue;


      // Compare with previous solutions: best or (distinct) second best solution?

      if (bestIndex == -1)
      { bestPValue = currentPValue;
        bestZScore = zScore;
        bestIndex = w;
      }
      else if ((COMPARISON_BY_PVALUE ? (zScore < bestZScore) : (properties.compareTo(listFeasibleAlignments.get(bestIndex)) < 0)))
      { if (! checkSolution(candidateLocation, candidateIsInvertedStrand, listFeasibleAlignments.get(bestIndex).location, listFeasibleAlignments.get(bestIndex).isInvertedStrand))
        { secondBestPValue = bestPValue;
          secondBestZScore = bestZScore;
          secondBestIndex = bestIndex;
        }
        else if (secondBestIndex != -1 && checkSolution(candidateLocation, candidateIsInvertedStrand, listFeasibleAlignments.get(secondBestIndex).location, listFeasibleAlignments.get(secondBestIndex).isInvertedStrand) )  // solution similar to the second one? delete the latter
	    { secondBestIndex = -1;
	      secondBestZScore = -4.0;
	      secondBestPValue = 1;
	    }
        bestPValue = currentPValue;
        bestZScore = zScore;
        bestIndex = w;
      }
      else if (secondBestIndex == -1 || (COMPARISON_BY_PVALUE ? (zScore < secondBestZScore) : (properties.compareTo(listFeasibleAlignments.get(secondBestIndex)) < 0)))
      { if (! checkSolution(candidateLocation, candidateIsInvertedStrand, listFeasibleAlignments.get(bestIndex).location, listFeasibleAlignments.get(bestIndex).isInvertedStrand))
        { secondBestPValue = currentPValue;
          secondBestZScore = zScore;
          secondBestIndex = w;
        }
        
      }
      

      
    }
    
    
        
    

    // Check for uniqueness:
    
    boolean isUnique = true;    

    boolean significantBestSolution = false;
    if ( bestIndex != -1 && ( COMPARISON_BY_PVALUE ? (bestPValue <= Fragment.PVALUE && listFeasibleAlignments.size() > Fragment.MINIMUMNUMBEROFFEASIBLESOLUTIONS) : (listFeasibleAlignments.get(bestIndex).score > RestrictionMap.MINIMUMSCORE) ) )

    { bestAlignment = listFeasibleAlignments.get(bestIndex);
      bestAlignment.secondBestPValue = secondBestPValue;
      bestAlignment.secondBestLocation = (secondBestPValue < 1 ? listFeasibleAlignments.get(secondBestIndex).location : -1);
       
       
      // Uniqueness test:
      
      if (secondBestIndex == -1)
      { significantBestSolution = true; }
      else
      { if ( COMPARISON_BY_PVALUE )
		{ if (secondBestPValue > Fragment.DIFFERENCESECONDBESTPVALUE * bestPValue)
		  { significantBestSolution = true;
          }
		  else
		  { isUnique = false;
		  }
		}
		else // F-test
		{ int bestLocation = bestAlignment.location;
	      boolean bestIsInvertedStrand = bestAlignment.isInvertedStrand;
	      
	      double bestScore = bestAlignment.score;
	      int bestMapMatches = bestAlignment.mapMatches - 2;
	      double bestWholeChiSquare = bestAlignment.chiSquare;
		
		
		  for (int w = 0; w < listFeasibleAlignments.size(); w++)
		  { if (w == bestIndex)
		      continue;
		   
		    CandidateSolution properties = listFeasibleAlignments.get(w);
		    int candidateLocation = properties.location;
		    boolean candidateIsInvertedStrand = properties.isInvertedStrand;
		
		    if (checkSolution(candidateLocation, candidateIsInvertedStrand, bestLocation, bestIsInvertedStrand))
		      continue;
		    
		    double candidateScore = properties.score;
		    int candidateMapMatches = properties.mapMatches - 2; // i.e., minus the end matches not counted towards the chi square test
		    double candidateWholeChiSquare = properties.chiSquare;
		
		    if (Math.abs(bestScore - candidateScore) < (2 * program2DP.matchBonus)) // => same number of cut errors, or at most 1 cut error difference
		    { double fval = betai((double)(bestMapMatches - 1) / 2.0, (double)(candidateMapMatches - 1) / 2.0, (bestWholeChiSquare + 2)/(bestWholeChiSquare + candidateWholeChiSquare + 4));
		      if ( errorOnMAXIT )
		      { System.err.println("errorOnMAXIT: " + opticalMap.index + " - " + opticalMap.id);
		        errorOnMAXIT = false;
		        isUnique = false;
		        break;
		      }
		       
		      if (fval >= FDISTRIBPVALUE)
		      { isUnique = false;
		        break;
		      }
		     }
		   }
			   
		   if (! isUnique)
			 significantBestSolution = true;
				
         }
      }
    }
                  	
    boolean notAligned = false;
 
    bestAlignment = null;

    if ((bestIndex != -1) &&  significantBestSolution )
    { bestAlignment = listFeasibleAlignments.get(bestIndex);
    }
    else
    { noAlignmentFound++;
      
      noNotAlignedMaps++;

      notAligned = true;
      
      if (bestIndex != -1)
      { bestAlignment = listFeasibleAlignments.get(bestIndex);
        bestAlignment.isUnique = isUnique;
        
        int seedsInSameLocation = 0;
        
        bestAlignment.seedsInSameLocation = seedsInSameLocation;

        bestAlignment.avgFragmentSizeOM = bestAlignment.avgFragmentSizeOM / (double)bestAlignment.countFragmentsOM;
        bestAlignment.avgFragmentSizeOMMinusEnds = bestAlignment.avgFragmentSizeOMMinusEnds / (double)bestAlignment.countFragmentsOMMinusEnds;
        bestAlignment.avgFragmentSizeReferenceMinusEnds = bestAlignment.avgFragmentSizeReferenceMinusEnds / (double)bestAlignment.countFragmentsReferenceMinusEnds;
        bestAlignment.fragSizeErrorAbs = bestAlignment.fragSizeErrorAbs / bestAlignment.mapMatches;
        bestAlignment.fragSizeErrorRelative = bestAlignment.fragSizeErrorRelative / bestAlignment.mapMatches;
        bestAlignment.secondBestPValue = secondBestPValue;
        bestAlignment.secondBestLocation = (secondBestPValue < 1 ?listFeasibleAlignments.get(secondBestIndex).location : -1);
        
        plot(opticalMap,bestAlignment, significantBestSolution,feasibleAlignments,seedsTried);
      }
      else
      { plotNotFoundSolution(opticalMap,false,feasibleAlignments,seedsTried);
      }
      
    }


    if (!notAligned)
    { // M: No. of solutions found similar to the best one (true or near-true solutions)
      
      int seedsInSameLocation = 0;

      bestAlignment.seedsInSameLocation = seedsInSameLocation;
      bestAlignment.isUnique = isUnique;

      bestAlignment.avgFragmentSizeOM = bestAlignment.avgFragmentSizeOM / (double)bestAlignment.countFragmentsOM;
      bestAlignment.avgFragmentSizeOMMinusEnds = bestAlignment.avgFragmentSizeOMMinusEnds / (double)bestAlignment.countFragmentsOMMinusEnds;
      bestAlignment.avgFragmentSizeReferenceMinusEnds = bestAlignment.avgFragmentSizeReferenceMinusEnds / (double)bestAlignment.countFragmentsReferenceMinusEnds;
      bestAlignment.fragSizeErrorAbs = bestAlignment.fragSizeErrorAbs / bestAlignment.mapMatches;
      bestAlignment.fragSizeErrorRelative = bestAlignment.fragSizeErrorRelative / bestAlignment.mapMatches;


      bestAlignment.significant = true;
      

      // Statistics

      plot(opticalMap,bestAlignment, significantBestSolution,feasibleAlignments,seedsTried);
      
    }

    if ( notAligned )
      return false; 

    noAlignedMaps++;
    
    return true;
    
  }



  int counter = 0;

  // M: loadReferenceSeeds  
 


  // Load optical maps

  public void loadOMs(String fileName, int firstIndex, int lastIndex) throws IOException
  { noOpticalMaps = 0;
    int fragmentSize;
    int location;
    int sumOfFragments = 0;
    discardedMaps = 0;
    int retainedMaps = 0;

    boolean allMaps = (firstIndex == -1 && lastIndex == -1);

    ArrayList<RestrictionMap> listOpticalMaps = (!allMaps ? new ArrayList<RestrictionMap>(lastIndex - firstIndex + 1) : new ArrayList<RestrictionMap>(100000));

    Scanner file = new Scanner(new File(fileName));
    Scanner line;

    while ( file.hasNextLine() ) // for each optical map
    { // Read heading:
      RestrictionMap opticalMap = new RestrictionMap();
      ArrayList<Fragment> opticalMapArrayList = new ArrayList<Fragment>(40);
      String id = file.nextLine();  // mapHeading

      // read fragment sizes:
      line = new Scanner(file.nextLine()).useDelimiter("\t");
      line.next(); //enzyme name
      line.next(); //'K'
      
      //if (noOpticalMaps == Fragment.MAXOMSTOLOAD) break;

      location = 0;
      sumOfFragments = 0;
      
      if ( (!allMaps && !(noOpticalMaps >= firstIndex && noOpticalMaps <= lastIndex)) )
      { noOpticalMaps++;
      }
      else
      { // Load fragment sizes of a single OM:
        while ( line.hasNext() )
        { fragmentSize = (int)((Float.parseFloat(line.next())) * 1000 + 0.5);
        
          opticalMapArrayList.add(new Fragment(fragmentSize, opticalMap, location));
  	      location++;
          sumOfFragments += fragmentSize;
        } // END WHILE
        
        
        // Skip map if its number of fragments or size are lower than a certain threshold:
        if (opticalMapArrayList.size() >= Fragment.LOWNUMBEROFFRAGMENTS && sumOfFragments >= Fragment.MINMOLECULESIZE)
        { noOpticalMaps++;
  
          opticalMap.fragments = new Fragment[location];
          opticalMapArrayList.toArray(opticalMap.fragments);
         
          opticalMap.noFragments = location;
          opticalMap.id = id;
          opticalMap.sizeMapInBases = sumOfFragments;
          opticalMap.index = noOpticalMaps; //Integer.parseInt(id.substring(4,id.length()));  // Format:  Map_XXX
  
          listOpticalMaps.add(opticalMap);
  
          // Account for the total size of the selected maps:
          totalSizeMaps += sumOfFragments;
          
          retainedMaps++;
        }
        else
        {
          
          noOpticalMaps++;
          
          discardedMaps++;
          
          opticalMap.fragments = new Fragment[location];
          opticalMapArrayList.toArray(opticalMap.fragments);
         
          opticalMap.noFragments = location;
          opticalMap.id = id;
          opticalMap.sizeMapInBases = sumOfFragments;
          opticalMap.index = noOpticalMaps;

          plotDiscardedMap(opticalMap);
        }
      }

      if ( file.hasNextLine() ) file.nextLine();  // empty row between optical maps
    }  //END WHILE FOR EACH OPTICAL MAP
    file.close();

    originalNoOpticalMaps = retainedMaps + discardedMaps;
    noOpticalMaps = retainedMaps;

    opticalMaps = new RestrictionMap[noOpticalMaps];
    listOpticalMaps.toArray(opticalMaps);
    
    System.err.println("Number of discarded maps: " + discardedMaps);
  }




  public boolean checkSolution (int location, boolean isInverted, int rightStartingLocation, boolean rightisInvertedStrand)
  { if (rightStartingLocation == -1)
      return false;
    if (Math.abs(location - rightStartingLocation) < 6 && isInverted == rightisInvertedStrand)
      return true;
    return false;
  }
 

  
  
  public void plot (RestrictionMap opticalMap, CandidateSolution bestAlignment, boolean significant, int noFeasibleAlignments, int seedsTried)
  {  
  
    double chiSquareSeedMatch = (bestAlignment.yOM.size - bestAlignment.zRef.size) / (double)bestAlignment.zRef.stDev;
    double chiSquareSeedMatchPow = Math.pow(chiSquareSeedMatch, 2.0);
    double relativeErrorSizeSeedMatch = (bestAlignment.yOM.size - bestAlignment.zRef.size) / (double)bestAlignment.zRef.size;
  
    String str = opticalMap.index + "\t" + opticalMap.id + "\t" + significant + "\t" + bestAlignment.isUnique + "\t" + opticalMap.size() + "\t" + opticalMap.sizeMapInBases + "\t" + bestAlignment.pValue + "\t" + bestAlignment.isInvertedStrand + "\t" + bestAlignment.location + "\t" + bestAlignment.lengthInReference + "\t" + bestAlignment.usefulSizeComparison + "\t" + bestAlignment.usefulOMSize + "\t" + bestAlignment.usefulRefSize + "\t" + (float)(bestAlignment.lengthInReference / (double)opticalMap.size()) + "\t" + (float)(bestAlignment.countFragmentsOMMinusEnds / (double)opticalMap.size()) + "\t" + bestAlignment.mapMatches + "\t" + bestAlignment.missingCuts + "\t" + bestAlignment.falseCuts + "\t" + bestAlignment.missingFragments + "\t" + (float)(bestAlignment.missingCuts / (double)(bestAlignment.lengthInReference - bestAlignment.missingFragments - bestAlignment.verySmallFragmentsNotDetectedAsMissing - bestAlignment.verySmallFragmentsDetectedAsMissing - 1.0) * 100.0) + "\t" + (float)(bestAlignment.falseCuts / (bestAlignment.refSizeMinusMissingFragments / 100000.0)) + "\t" + (bestAlignment.noSmallRefFragments > 0 ? (float)((bestAlignment.missingFragments + bestAlignment.verySmallFragmentsDetectedAsMissing) / (double)bestAlignment.noSmallRefFragments  * 100.0) : 0) + "\t" + bestAlignment.chiSquarePerMatch + "\t" + (bestAlignment.chiSquare / bestAlignment.mapMatches) + "\t" + bestAlignment.chiSquare + "\t" + bestAlignment.score + "\t" + bestAlignment.probMissingCuts + "\t" + bestAlignment.probFalseCuts + "\t" + bestAlignment.probChiSquarePerMatch + "\t" + bestAlignment.zScore + "\t" + bestAlignment.zScoreEC + "\t" + bestAlignment.zScoreFCMM + "\t" + bestAlignment.zScoreCSPM + "\t" + bestAlignment.avgFragmentSizeOM + "\t" + bestAlignment.avgFragmentSizeOMMinusEnds + "\t" + bestAlignment.avgFragmentSizeReferenceMinusEnds + "\t" + bestAlignment.fragSizeErrorAbs + "\t" + bestAlignment.fragSizeErrorRelative + "\t"     + bestAlignment.countFragmentsOM + "\t" + bestAlignment.countFragmentsOMMinusEnds + "\t" + bestAlignment.countFragmentsReference + "\t" + bestAlignment.countFragmentsReferenceMinusEnds + "\t" + bestAlignment.noSmallRefFragments + "\t" + bestAlignment.refSizeMinusMissingFragments + "\t" + bestAlignment.yOM.location + "\t" + bestAlignment.yOM.size + "\t" + bestAlignment.zRef.size + "\t" + (Math.abs(chiSquareSeedMatch)) + "\t" + chiSquareSeedMatchPow + "\t" + relativeErrorSizeSeedMatch + "\t" + bestAlignment.seedsInSameLocation + "\t" +  noFeasibleAlignments + "\t" + seedsTried + "\t" + bestAlignment.verySmallFragmentsNotDetectedAsMissing + "\t" + bestAlignment.verySmallFragmentsDetectedAsMissing + "\t" + bestAlignment.secondBestPValue + "\t" + bestAlignment.secondBestLocation;
  
     if ( significant )
       outOK.println(str);
     else
       outNotOK.println(str);
  }
  
  public void plotNotFoundSolution (RestrictionMap opticalMap, boolean significant, int noFeasibleAlignments, int seedsTried)
  { String str = opticalMap.index + "\t" + opticalMap.id + "\t" + significant + "\t" + opticalMap.size() + "\t" + opticalMap.sizeMapInBases + "\t" + noFeasibleAlignments + "\t" + seedsTried;  
      
    outNotFound.println(str);
  }

  
  public void plotDiscardedMap (RestrictionMap opticalMap)
  { String str = opticalMap.index + "\t" + opticalMap.id + "\t" + opticalMap.size() + "\t" + opticalMap.sizeMapInBases;
      
    outDiscarded.println(str);
  }
  

  public double gammln(double xx)                // static
  {
    double x,y,tmp,ser;
    double[] cof = {76.18009172947146,-86.50532032941677,
 			24.01409824083091,-1.231739572450155,
  			0.1208650973866179e-2,-0.5395239384953e-5};
    int j;
    y = xx;
    x = xx;
    tmp= x+5.5;
    tmp -= (x+0.5)*Math.log(tmp);
    ser = 1.000000000190015;
    for (j = 0; j <= 5; j++)
    { y++;
      ser += cof[j] / y;
    }
    return (Math.log(2.5066282746310005 * ser / x) - tmp);
  }


  public double betacf (double a, double b, double x)     // static
  {
    int MAXIT = 100;
    double EPS = 3.0e-7;
    double FPMIN = 1.0e-30;
  
    int m,m2;
    double aa,c,d,del,h,qab,qam,qap;
    qab=a+b;
    qap=a+1.0;
    qam=a-1.0;
    c=1.0;
    d= 1.0 - (qab * x / qap);
  
    if (Math.abs(d) < FPMIN) d=FPMIN;
  
    d = 1.0 / d;
    h=d;
    for (m=1; m <= MAXIT; m++) {
  
      m2=2*m;
      aa=m*(b-m)*x/((qam+m2)*(a+m2));
      d=1.0+aa*d;
      
      if (Math.abs(d) < FPMIN) d=FPMIN;
      c=1.0+aa/c;
      if (Math.abs(c) < FPMIN) c=FPMIN;
  
      d=1.0/d;
      h *= (d * c);
      aa = (-(a+m))*(qab+m)*x/((a+m2)*(qap+m2));
      d=1.0+aa*d; 
  
      if (Math.abs(d) < FPMIN) d=FPMIN;
      c=1.0+aa/c;
      if (Math.abs(c) < FPMIN) c=FPMIN;
  
      d = 1.0 / d;
      del = d * c;
      h *= del;
  
      if (Math.abs(del - 1.0) < EPS) break;
    }
  
    if (m > MAXIT)
    { System.err.println("ERROR:  a or b too big, or MAXIT too small in routine 'betacf'");
      System.err.println(a + " " + b + " " + x);
      errorOnMAXIT = true;
    }
  
    return h;
  }
  

  public double betai (double a, double b, double x)     // static
  {
    double bt;
    if (x < 0.0 || x > 1.0)
    { System.err.println("ERROR:  Bad x in routine 'betai'");
    } //return -1;

    if (x == 0.0 || x == 1.0)
    { bt = 0.0; }
    else
    { bt = Math.pow(Math.E, this.gammln(a + b) - this.gammln(a) - this.gammln(b)+ a*Math.log(x) + b*Math.log(1.0 - x)); }
    
    if (x < ((a + 1.0)/(a + b + 2.0)))
    { return (bt * this.betacf(a, b, x) / a); }

    return (1.0 - bt * this.betacf(b, a, 1.0 - x) / b);
  }

}
