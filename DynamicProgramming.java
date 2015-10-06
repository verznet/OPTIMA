/*  DynamicProgramming v4

    Dynamic programming algorithm for aligning an entire optical map against
    part of a reference map, starting from a specific location or seed.
    
    
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
import java.util.ArrayList;
import java.util.Collections;

import java.io.Writer;
import java.io.PrintWriter;
import java.io.FileWriter;

import org.apache.commons.math3.distribution.*;
import org.apache.commons.math3.random.EmpiricalDistribution;
import org.apache.commons.math3.exception.NotStrictlyPositiveException;


public class DynamicProgramming
{ 
  public class Pair
  {  double bestScore;
     int bestLastLocationJ;

     public Pair(double aScore, int aLocation)
     {  bestScore = aScore;
        bestLastLocationJ = aLocation;
     }
  }



  int TOTALFRAGMENTS = 450000;
  Fragment[] reference = new Fragment[TOTALFRAGMENTS];

  int noReferenceFragments = 0;

  int maxNoOfConsecutiveFragmentsAlignedTogetherI;
  int maxNoOfConsecutiveFragmentsAlignedTogetherJ;
  double MAXTOTALJOINEDFRAGMENTSALLOWED;
  int matchBonus = (int)(RestrictionMap.maxNoOfMatches * Math.pow(Fragment.CSIGMA, 2.0));
  int matchBonus2 = 250 * matchBonus;
  int mapMatches = 0;
  long sumOfRefFragments = 0;
  int[][] matches;
  int[][] directions;
  int[][] directionsI;
  int[][] directionsJ;
  double maxChiSquare = 0;
  int locationMaxScoreJ = 0;
  boolean SOMAV3 = false;
  
  int avgFragmentSizeOM = 0;
  int avgFragmentSizeOMMinusEnds = 0;
  int avgFragmentSizeReferenceMinusEnds = 0;
  double fragSizeErrorAbs = 0.0;
  double fragSizeErrorRelative = 0.0;
  
  int countFragmentsOM = 0;
  int countFragmentsOMMinusEnds = 0;
  int countFragmentsReference = 0;
  int countFragmentsReferenceMinusEnds = 0;

  int totalJoinedFragmentsPercentageRef = 0;
  int totalJoinedFragmentsPercentageOM = 0;
  
  int currentMissingCuts = 0;
  int currentFalseCuts = 0;
  int currentMissingFragments = 0;
  
  int verySmallFragmentsNotDetectedAsMissing = 0;
  int verySmallFragmentsDetectedAsMissing = 0;

  int currentNoSmallRefFragments = 0; 
  int currentRefSizeMinusMissingFragments = 0;
  int currentUsefulRefSize = 0;
  int currentUsefulOMSize = 0;

  Fragment[] listReferenceSeedsToArray;
  ArrayList<Fragment> listReferenceSeeds;
  int noReferenceSeeds = 0;
  
  
  // *** Used in Fragment sizing error estimation and Violin plots:
  double[][] localChiSquare;
  int[][] localOMSize;
  int[][] localRefSize;
    

  
  public Pair dynamicProgrammingSOMAv2 (Fragment[] opticalMap, int locationInReference, boolean reverse, int highestLocationInReference)
  { // given a seed, align the map to the reference maps by extending the seed
    
    if (opticalMap.length == 0) return new Pair(0, -1);

    int lowerLevel = locationInReference;
    int upperBound = 0;

    if (opticalMap.length < 3) //4)
      upperBound = (int)(opticalMap.length * Fragment.BOUNDARYBIG + Fragment.BOUNDARYBIGCONSTANT);
    else
      upperBound = (int)(opticalMap.length * Fragment.BOUNDARYSMALL + Fragment.BOUNDARYSMALLCONSTANT);

    // Check boundaries:
    if (! reverse)
    { if (lowerLevel + upperBound - 1 > highestLocationInReference)
      { upperBound = highestLocationInReference - lowerLevel + 1; }
    }
    else
    { if ((lowerLevel - upperBound + 1) < 0)
      { upperBound = lowerLevel + 1; }
    }

    double[][] scores = new double[opticalMap.length + 1][upperBound + 1];

    int[] sumSizesI = new int[opticalMap.length + 1];
    long[] sumStDevJ = new long[upperBound + 1];
    int[] sumSizesJ = new int[upperBound + 1];
    
    int[][] matches = new int[opticalMap.length + 1][upperBound + 1];
    int[][] missingFragments = new int[opticalMap.length + 1][upperBound + 1];
    localChiSquare = new double[opticalMap.length + 1][upperBound + 1];
    localOMSize = new int[opticalMap.length + 1][upperBound + 1];
    localRefSize = new int[opticalMap.length + 1][upperBound + 1];
    int[][] localRefSizeWithMissing = new int[opticalMap.length + 1][upperBound + 1];
    int[][] refSizeMinusMissingFragments = new int[opticalMap.length + 1][upperBound + 1];

    directionsI = new int[opticalMap.length + 1][upperBound + 1];
    directionsJ = new int[opticalMap.length + 1][upperBound + 1];

    int i;
    int j;
    int k;
    int l;
    double maxGlobalScore = RestrictionMap.MINIMUMSCORE;
    int locationMaxScoreI = 0;
    locationMaxScoreJ = 0;

    // Initialization:
    sumSizesJ[0] = 0;
    sumStDevJ[0] = 0;
    
    for (j = 1; j <= upperBound; j++)
    { if ( reverse )
      { sumSizesJ[j] = sumSizesJ[j - 1] + reference[lowerLevel - j + 1].size;
        sumStDevJ[j] = sumStDevJ[j - 1] + reference[lowerLevel - j + 1].squaredStDev;
      }
      else
      { sumSizesJ[j] = sumSizesJ[j - 1] + reference[lowerLevel + j - 1].size;
        sumStDevJ[j] = sumStDevJ[j - 1] + reference[lowerLevel + j - 1].squaredStDev;
      }
    }

    sumSizesI[0] = 0;

    for (i = 1; i <= opticalMap.length; i++)
    { sumSizesI[i] = sumSizesI[i - 1] + opticalMap[i - 1].size;
    }

    scores[0][0] = 0.01;
    matches[0][0] = 0;
    directionsI[0][0] = 0;
    directionsJ[0][0] = 0;

    double chiSquare = 0;

    int POSSIBLEFEASIBLESOLUTIONS = Fragment.maxNoOfConsecutiveFragmentsAlignedTogetherIFlexible * upperBound + 3;
    int numberOfNoFeasibleSolutionSteps = 0;

    boolean thereIsAScore = false;

MAINFOR: for (i = 1; i <= opticalMap.length; i++)
    { int jLeftBound = i - 12;
      if (jLeftBound < 1) jLeftBound = 1;

      int jRightBound = upperBound;

      for (j = jLeftBound; j <= jRightBound; j++)
      { double currentScore = RestrictionMap.MINIMUMSCORE;
        double maxScore = RestrictionMap.MINIMUMSCORE;

        directionsI[i][j] = -1;
        directionsJ[i][j] = -1;

        for (k = i; (k >= 1) && (k > (i - maxNoOfConsecutiveFragmentsAlignedTogetherI)); k--)   // 1...i
        {
          for (l = j; (l >= 1) && (l > (j - maxNoOfConsecutiveFragmentsAlignedTogetherJ)); l--)     // 1...j
          {
            double previousScore = scores[k-1][l-1];
            currentScore = 0;
            boolean gotBestSolutionPreviouslyWithLessGaps = false;

            if (previousScore <= RestrictionMap.MINIMUMSCORE || (previousScore == 0))
            { continue; }

            // Limit DP to match the entire OM with the initial Reference fragments
            if (!(((l - 1) >= 1 && (k - 1) >= 1)  ||  ((l - 1) == 0 && (k - 1) == 0)))
            { continue; }

            int otScore = sumSizesI[i] - sumSizesI[k - 1];

            int csScore = sumSizesJ[j] - sumSizesJ[l - 1];
            int sigmasScore = (int)(sumStDevJ[j] - sumStDevJ[l - 1]);
            
            int sizeGapJ = sumSizesJ[l] - sumSizesJ[l - 1];      // gap on l
            int sizeGapJsecond = (l+1 < j ? sumSizesJ[l+1] - sumSizesJ[l] : 0);
            int sizeGapJthird = (l+2 < j ? sumSizesJ[l+2] - sumSizesJ[l+1] : 0);

            // Try with 1 gap at the beginning of the match
            int csScoreWithGap = sumSizesJ[j] - sumSizesJ[l];
            int sigmasScoreWithGap = (int)(sumStDevJ[j] - sumStDevJ[l]);

            double chiSquare2 = (otScore - csScore) / Math.sqrt((double)sigmasScore);
            chiSquare = Math.abs(chiSquare2);
            double chiSquarePow = Math.pow(chiSquare2, 2.0);

            double chiSquareWithGapJ2 = (otScore - csScoreWithGap) / Math.sqrt((double)sigmasScoreWithGap);
            double chiSquareWithGapJ = Math.abs(chiSquareWithGapJ2);
            double chiSquareWithGapJPow = Math.pow(chiSquareWithGapJ2, 2.0);

            int currentMatches = matches[k - 1][l - 1] + 1;

            double missingCutsPercentageRef = (j - currentMatches);
            double falseCutsPercentageOM = (i - currentMatches);
            double nonMissingFragmentsCuts = j - missingFragments[k-1][l-1] - 1.0; // number of cuts minus missing fragments                    


            if((i == opticalMap.length && chiSquare2 <= Fragment.CSIGMA) || chiSquare <= Fragment.CSIGMA)
            { if ( SOMAV3 )
              { currentScore = (-matchBonus2 * (i-k)) - (matchBonus * (j-l)) - chiSquarePow + previousScore;
              }
              else
              { if (i == opticalMap.length) currentScore = -matchBonus * (j-l + Fragment.FALSECUTSSCORE*(i-k)) + previousScore;
                else currentScore = -matchBonus * (j-l + Fragment.FALSECUTSSCORE*(i-k)) - chiSquarePow + previousScore;
              }

              if ( maxScore <= RestrictionMap.MINIMUMSCORE || (currentScore > maxScore  &&  (localChiSquare[i][j] == 0.0 || chiSquarePow < Fragment.CHISQUAREREWARDCONSTANT * localChiSquare[i][j])) || (localChiSquare[i][j] != 0.0 && Math.abs(currentScore - maxScore) < Fragment.SIMILARITYSCORES && localChiSquare[i][j] > Fragment.CHISQUAREREWARDCONSTANT * chiSquarePow) )              { 
                // M: Check first:
                
                maxScore = currentScore;
                matches[i][j] = currentMatches;
                directionsI[i][j] = k - 1;
                directionsJ[i][j] = l - 1;
                missingFragments[i][j] = missingFragments[k-1][l-1];
                if (i == opticalMap.length)
                { localChiSquare[i][j] = 0.0;
                  localOMSize[i][j] = 0;
                  localRefSize[i][j] = 0;
                  localRefSizeWithMissing[i][j] = 0;
                  
                  refSizeMinusMissingFragments[i][j] = (int)(otScore / Fragment.ADJUSTMENT_OPTICAL_MAP_SIZE) + 1;    // !!!
                }
                else
                { localChiSquare[i][j] = chiSquarePow;
                  localOMSize[i][j] = otScore;
                  localRefSize[i][j] = csScore;
                  localRefSizeWithMissing[i][j] = csScore;
                  
                  refSizeMinusMissingFragments[i][j] = csScore;
                }

                gotBestSolutionPreviouslyWithLessGaps = true;
              }

            }
            

            //else allow 1 initial mismatch within the match
            if (l < j && sizeGapJ < Fragment.MAXSIZESMALLFRAGMENTS && ( (i == opticalMap.length && chiSquareWithGapJ2 <= Fragment.CSIGMA) || (chiSquareWithGapJ <= Fragment.CSIGMA) ) ) // || chiSquareWithGapI <= Fragment.CSIGMA))                     & k < i      i < opticalMap.length &&
            { if ( SOMAV3 )
              { currentScore = (-matchBonus2 * (i-k)) - (matchBonus * (j-l)) - chiSquareWithGapJPow + previousScore;
              }
              else
              { if (i == opticalMap.length) currentScore = -matchBonus * (j-l + Fragment.FALSECUTSSCORE*(i-k)) + previousScore;
                else currentScore = -matchBonus * (j-l + Fragment.FALSECUTSSCORE*(i-k)) - chiSquareWithGapJPow + previousScore;
              }

              if (maxScore <= RestrictionMap.MINIMUMSCORE   ||   (currentScore > maxScore  &&  (localChiSquare[i][j] == 0.0 || chiSquareWithGapJPow < Fragment.CHISQUAREREWARDCONSTANT * localChiSquare[i][j])) || (localChiSquare[i][j] != 0.0 && Math.abs(currentScore - maxScore) < Fragment.SIMILARITYSCORES && localChiSquare[i][j] > Fragment.CHISQUAREREWARDCONSTANT * chiSquareWithGapJPow) || (gotBestSolutionPreviouslyWithLessGaps  && localChiSquare[i][j] != 0.0 && Math.abs(currentScore - maxScore) < Fragment.SIMILARITYSCORES && localChiSquare[i][j] > Fragment.CHISQUAREREWARDCONSTANTFORGAPS * chiSquareWithGapJPow ) )
              { 
                // M: Check first:
                double currentGap = 1.0;
                
                maxScore = currentScore;
                matches[i][j] = currentMatches;
		        directionsI[i][j] = k - 1;
                directionsJ[i][j] = l - 1;
                missingFragments[i][j] = missingFragments[k-1][l-1] + 1;
                if (i == opticalMap.length)
                { localChiSquare[i][j] = 0.0;
                  localOMSize[i][j] = 0;
                  localRefSize[i][j] = 0;
                  localRefSizeWithMissing[i][j] = 0;
                  
                  refSizeMinusMissingFragments[i][j] = (int)(otScore / Fragment.ADJUSTMENT_OPTICAL_MAP_SIZE) + 1;    // !!!
                }
                else
                { localChiSquare[i][j] = chiSquareWithGapJPow;
                  localOMSize[i][j] = otScore;
                  localRefSize[i][j] = csScoreWithGap;
                  localRefSizeWithMissing[i][j] = csScore;

                  refSizeMinusMissingFragments[i][j] = csScoreWithGap;
                }

                gotBestSolutionPreviouslyWithLessGaps = true;
              }
            }
            
            //else allow 2 initial mismatches within the match
            if (l+1 < j && sizeGapJ < Fragment.MAXSIZESMALLFRAGMENTS && sizeGapJsecond < Fragment.MAXSIZESMALLFRAGMENTS)
            { int csScoreWithGap2 = sumSizesJ[j] - sumSizesJ[l+1];
              int sigmasScoreWithGap2 = (int)(sumStDevJ[j] - sumStDevJ[l+1]);

              double chiSquareWithGapJJ2 = (otScore - csScoreWithGap2) / Math.sqrt((double)sigmasScoreWithGap2);
              double chiSquareWithGapJJ = Math.abs(chiSquareWithGapJJ2);
              double chiSquareWithGapJJPow = Math.pow(chiSquareWithGapJJ2, 2.0);

              if ( !((i == opticalMap.length && chiSquareWithGapJJ2 <= Fragment.CSIGMA) || (chiSquareWithGapJJ <= Fragment.CSIGMA)) )
              { }
              else
              {
                if ( SOMAV3 )
                { currentScore = (-matchBonus2 * (i-k)) - (matchBonus * (j-l)) - chiSquareWithGapJJPow + previousScore;
                }
                else
                { if (i == opticalMap.length) currentScore = -matchBonus * (j-l + Fragment.FALSECUTSSCORE*(i-k)) + previousScore;
                  else currentScore = -matchBonus * (j-l + Fragment.FALSECUTSSCORE*(i-k)) - chiSquareWithGapJJPow + previousScore;
                }
  
                if (maxScore <= RestrictionMap.MINIMUMSCORE || (currentScore > maxScore  &&  (localChiSquare[i][j] == 0.0 || chiSquareWithGapJJPow < Fragment.CHISQUAREREWARDCONSTANTFORGAPS * localChiSquare[i][j])) || (localChiSquare[i][j] != 0.0 && Math.abs(currentScore - maxScore) < Fragment.SIMILARITYSCORES && localChiSquare[i][j] > Fragment.CHISQUAREREWARDCONSTANTFORGAPS * chiSquareWithGapJJPow) ||( gotBestSolutionPreviouslyWithLessGaps  && localChiSquare[i][j] != 0.0 && Math.abs(currentScore - maxScore) < Fragment.SIMILARITYSCORES && localChiSquare[i][j] > Fragment.CHISQUAREREWARDCONSTANTFORGAPS * chiSquareWithGapJJPow))
                { 
                  // M: Check first:
                  double currentGap = 2.0;
                    
                  maxScore = currentScore;
                  matches[i][j] = currentMatches;
  		directionsI[i][j] = k - 1;
                  directionsJ[i][j] = l - 1;
                  missingFragments[i][j] = missingFragments[k-1][l-1] + 2;
                  if (i == opticalMap.length)
                  { localChiSquare[i][j] = 0.0;
                    localOMSize[i][j] = 0;
                    localRefSize[i][j] = 0;
                    localRefSizeWithMissing[i][j] = 0;
                    
                    refSizeMinusMissingFragments[i][j] = (int)(otScore / Fragment.ADJUSTMENT_OPTICAL_MAP_SIZE) + 1;    // !!!
                  }
                  else
                  { localChiSquare[i][j] = chiSquareWithGapJJPow;
                    localOMSize[i][j] = otScore;
                    localRefSize[i][j] = csScoreWithGap2;
                    localRefSizeWithMissing[i][j] = csScore;
                    
                    refSizeMinusMissingFragments[i][j] = csScoreWithGap2;
                  }
                                   
                  gotBestSolutionPreviouslyWithLessGaps = true;
                }
              }
            }


            //else allow 3 initial mismatches within the match
            if (l+2 < j && sizeGapJ < Fragment.MAXSIZESMALLFRAGMENTS && sizeGapJsecond < Fragment.MAXSIZESMALLFRAGMENTS && sizeGapJthird < Fragment.MAXSIZESMALLFRAGMENTS)
            { int csScoreWithGap3 = sumSizesJ[j] - sumSizesJ[l+2];
              int sigmasScoreWithGap3 = (int)(sumStDevJ[j] - sumStDevJ[l+2]);

              double chiSquareWithGapJJJ2 = (otScore - csScoreWithGap3) / Math.sqrt((double)sigmasScoreWithGap3);
              double chiSquareWithGapJJJ = Math.abs(chiSquareWithGapJJJ2);
              double chiSquareWithGapJJJPow = Math.pow(chiSquareWithGapJJJ2, 2.0);

              if ( !((i == opticalMap.length && chiSquareWithGapJJJ2 <= Fragment.CSIGMA) || (chiSquareWithGapJJJ <= Fragment.CSIGMA)) )
              { }
              else
              {
                if ( SOMAV3 )
                { currentScore = (-matchBonus2 * (i-k)) - (matchBonus * (j-l)) - chiSquareWithGapJJJPow + previousScore;
                }
                else
                { if (i == opticalMap.length) currentScore = -matchBonus * (j-l + Fragment.FALSECUTSSCORE*(i-k)) + previousScore;
                  else currentScore = -matchBonus * (j-l + Fragment.FALSECUTSSCORE*(i-k)) - chiSquareWithGapJJJPow + previousScore;
                }
  
                if (maxScore <= RestrictionMap.MINIMUMSCORE || (currentScore > maxScore  &&  (localChiSquare[i][j] == 0.0 || chiSquareWithGapJJJPow < Fragment.CHISQUAREREWARDCONSTANTFORGAPS * localChiSquare[i][j])) || (localChiSquare[i][j] != 0.0 && Math.abs(currentScore - maxScore) < Fragment.SIMILARITYSCORES && localChiSquare[i][j] > Fragment.CHISQUAREREWARDCONSTANTFORGAPS * chiSquareWithGapJJJPow) ||( gotBestSolutionPreviouslyWithLessGaps  && localChiSquare[i][j] != 0.0 && Math.abs(currentScore - maxScore) < Fragment.SIMILARITYSCORES && localChiSquare[i][j] > Fragment.CHISQUAREREWARDCONSTANTFORGAPS * chiSquareWithGapJJJPow))
                { 
                  // M: Check first:
                  double currentGap = 3.0;
                    
                  maxScore = currentScore;
                  matches[i][j] = currentMatches;
  		directionsI[i][j] = k - 1;
                  directionsJ[i][j] = l - 1;
                  missingFragments[i][j] = missingFragments[k-1][l-1] + 3;
                  if (i == opticalMap.length)
                  { localChiSquare[i][j] = 0.0;
                    localOMSize[i][j] = 0;
                    localRefSize[i][j] = 0;
                    localRefSizeWithMissing[i][j] = 0;
                    
                    refSizeMinusMissingFragments[i][j] = (int)(otScore / Fragment.ADJUSTMENT_OPTICAL_MAP_SIZE) + 1;    // !!!
                  }
                  else
                  { localChiSquare[i][j] = chiSquareWithGapJJJPow;
                    localOMSize[i][j] = otScore;
                    localRefSize[i][j] = csScoreWithGap3;
                    localRefSizeWithMissing[i][j] = csScore;
                    
                    refSizeMinusMissingFragments[i][j] = csScoreWithGap3;
                  }

                  gotBestSolutionPreviouslyWithLessGaps = true;
                }
              }
            }


          } //END FOR j
        } //END FOR k

        // Check if there are feasible solutions in the previous cells:
        scores[i][j] = maxScore;

        if (maxScore > RestrictionMap.MINIMUMSCORE)
        { numberOfNoFeasibleSolutionSteps = 0;
          thereIsAScore = true;
        }
        else 
        { numberOfNoFeasibleSolutionSteps++;
          if (numberOfNoFeasibleSolutionSteps > POSSIBLEFEASIBLESOLUTIONS)
          { break MAINFOR;
          }
        }
      }
    }
    
    if (! thereIsAScore)
    {
      return new Pair(RestrictionMap.MINIMUMSCORE, -1);
    }



    // Compute location max scores
    
    i = opticalMap.length;
    int lowerBound = opticalMap.length - Fragment.maxNoOfConsecutiveFragmentsAlignedTogetherIFlexible - 1;
    if (lowerBound < 1) lowerBound = 1;

Mx: for (j = lowerBound; j <= upperBound; j++)
    { if (scores[i][j] == 0 || scores[i][j] <= maxGlobalScore)
      { continue; }

      maxGlobalScore = scores[i][j];
      locationMaxScoreJ = j;
    }
    
    if (maxGlobalScore <= RestrictionMap.MINIMUMSCORE)
      return new Pair(RestrictionMap.MINIMUMSCORE, -1);
       
    j = locationMaxScoreJ;

    this.totalJoinedFragmentsPercentageRef = j - matches[i][j];
    this.totalJoinedFragmentsPercentageOM = i - matches[i][j];
    
    this.currentMissingFragments = missingFragments[i][j];
    this.currentMissingCuts = this.totalJoinedFragmentsPercentageRef - this.currentMissingFragments;
    this.currentFalseCuts = this.totalJoinedFragmentsPercentageOM;
    
    

    // Find back the optimal path to compute the statistical significance of matches
    
    i = opticalMap.length;
    j = locationMaxScoreJ;
    
    int previousI;
    int previousJ;
    
        
    // Compute error values      
    avgFragmentSizeOM = 0;
    avgFragmentSizeOMMinusEnds = 0;
    avgFragmentSizeReferenceMinusEnds = 0;
    fragSizeErrorAbs = 0.0;
    fragSizeErrorRelative = 0.0;
    
    countFragmentsOM = 0;
    countFragmentsOMMinusEnds = 0;
    countFragmentsReference = 0;
    countFragmentsReferenceMinusEnds = 0;
    

    currentNoSmallRefFragments = 0; 
    currentRefSizeMinusMissingFragments = 0;
    currentUsefulRefSize = 0;
    currentUsefulOMSize = 0;
    
    verySmallFragmentsNotDetectedAsMissing = 0; 


    for (int pos = 0; pos < locationMaxScoreJ; pos++)
    { if ( reverse )
      { if (reference[lowerLevel - pos].size < Fragment.MAXSIZESMALLFRAGMENTS)
          currentNoSmallRefFragments++;
      }
      else
      { if (reference[lowerLevel + pos].size < Fragment.MAXSIZESMALLFRAGMENTS)
          currentNoSmallRefFragments++;
      }
    }

    i = opticalMap.length;
    j = locationMaxScoreJ;
    int tmpSum = 0;
    
    int verySmallFragmentsDetectedAsMissing = 0;
    
    
    while (i > 0 && j > 0)
    { previousI = i;
      previousJ = j;

      if (localOMSize[i][j] != 0 && localRefSize[i][j] != 0)
      { currentUsefulRefSize += localRefSize[i][j];
        currentUsefulOMSize += localOMSize[i][j];
        
        avgFragmentSizeOMMinusEnds += localOMSize[i][j];
        avgFragmentSizeReferenceMinusEnds += localRefSizeWithMissing[i][j];
        
        fragSizeErrorAbs += (Math.abs(localOMSize[i][j] - localRefSize[i][j]) / (double)localRefSize[i][j]);
        fragSizeErrorRelative += ((localOMSize[i][j] - localRefSize[i][j]) / (double)localRefSize[i][j]);
      }

      avgFragmentSizeOM += localOMSize[i][j];

      currentRefSizeMinusMissingFragments += refSizeMinusMissingFragments[i][j];


      i = directionsI[previousI][previousJ];
      j = directionsJ[previousI][previousJ];

      countFragmentsOM += previousI - i;
      countFragmentsReference += previousJ - j;

      if (localOMSize[i][j] != 0 && localRefSize[i][j] != 0)
      { countFragmentsOMMinusEnds += previousI - i;
        countFragmentsReferenceMinusEnds += previousJ - j;
      }
      
      int noMissingFragmentsDetectedAtCurrentMatch = missingFragments[previousI][previousJ] - missingFragments[i][j]; // actually the previous match...
      
      int currentVerySmallFragmentsNotDetectedAsMissing = 0;

      
      if (previousJ > j + 1)
      { XFor: for (int prevBackward = j+1; prevBackward <= previousJ; prevBackward++)
        { if (noMissingFragmentsDetectedAtCurrentMatch > 0)
          { // first missing fragments detected for current match (up to 3 missing frags)
            noMissingFragmentsDetectedAtCurrentMatch--;
            
            int tmpSize = (reverse ? reference[lowerLevel - prevBackward + 1].size : reference[lowerLevel + prevBackward - 1].size);
            if (tmpSize <= Fragment.MAXSIZEVERYSMALLFRAGMENTS) 
            { verySmallFragmentsDetectedAsMissing++;
            }
          }
          else 
          { int tmpSize = (reverse ? reference[lowerLevel - prevBackward + 1].size : reference[lowerLevel + prevBackward - 1].size);
            if (tmpSize <= Fragment.MAXSIZEVERYSMALLFRAGMENTS)
            { currentVerySmallFragmentsNotDetectedAsMissing++;
              
            }
          }
        }
        
        noMissingFragmentsDetectedAtCurrentMatch = missingFragments[previousI][previousJ] - missingFragments[i][j];
        
        if (currentVerySmallFragmentsNotDetectedAsMissing == previousJ - j - noMissingFragmentsDetectedAtCurrentMatch)
          currentVerySmallFragmentsNotDetectedAsMissing--; // we expect to see at least a "non-small" fragment in a match!
        if (currentVerySmallFragmentsNotDetectedAsMissing < 0)
        { currentVerySmallFragmentsNotDetectedAsMissing = 0;
          System.err.println("Error on negative very small frags not detected as missing!");
        }
      }
      
      verySmallFragmentsNotDetectedAsMissing += currentVerySmallFragmentsNotDetectedAsMissing;
      
    }
    
    this.currentMissingFragments -= verySmallFragmentsDetectedAsMissing;
    this.currentMissingCuts = this.currentMissingCuts - verySmallFragmentsNotDetectedAsMissing;
    if (this.currentMissingCuts < 0) System.err.println("$$ MC:" + this.currentMissingCuts + " SmallMFnD:" + verySmallFragmentsNotDetectedAsMissing + " TotMF:" + this.currentMissingFragments + " TotJF:" + this.totalJoinedFragmentsPercentageRef + " LengthRef:" + locationMaxScoreJ);
    
    this.verySmallFragmentsDetectedAsMissing = verySmallFragmentsDetectedAsMissing;

    return new Pair(maxGlobalScore, ( reverse ? lowerLevel - (locationMaxScoreJ - 1) : lowerLevel + (locationMaxScoreJ - 1)));
  }



  public void loadReferenceFragments(String inSilicoMaps) throws IOException
  { long totalSize = 0;
    listReferenceSeeds = new ArrayList<Fragment>(TOTALFRAGMENTS*2);
    noReferenceFragments = 0;
    noReferenceSeeds = 0;
    int indexChromosome = 1;

    Scanner file = new Scanner(new File(inSilicoMaps));
    Scanner line = new Scanner(file.nextLine());
    String restrSite = line.next();
    
    while (file.hasNextLine() && restrSite.charAt(0) == '>') 
    {  int stepIndex = 0;     // change of chromosome

       int chrSize = 0;
       int restrictionSitePrev = 0;
       int restrictionSiteCurrent;
       int fragmentSize;
       
       // heading: read chr size
       //line.next(); // chr name  // already done
       chrSize = Integer.parseInt(line.next());
       line.next(); // number of restriction sites

       restrSite = "";

       while (file.hasNextLine() && (restrSite = (line = new Scanner(file.nextLine())).next()).charAt(0) != '>')
       { // restriction sites
      
         do
         { restrictionSiteCurrent = Integer.parseInt(restrSite);
           fragmentSize = restrictionSiteCurrent - restrictionSitePrev;
      
           reference[noReferenceFragments] = new Fragment(fragmentSize, noReferenceFragments, totalSize, false);

           Fragment fragNorm = new Fragment(fragmentSize, noReferenceFragments, totalSize, false);
           Fragment fragRev = new Fragment(fragmentSize, noReferenceFragments, totalSize, true);
           listReferenceSeeds.add(fragNorm);
           listReferenceSeeds.add(fragRev);
           noReferenceSeeds += 2;
           noReferenceFragments++;
           stepIndex++;
          	
           totalSize =+ fragmentSize;

           restrictionSitePrev = restrictionSiteCurrent;
         } while (line.hasNext() && (restrSite = line.next()) != null);   // END DO WHILE
       } // END WHILE
      
       // Last fragment:
       restrictionSiteCurrent = chrSize;
       fragmentSize = restrictionSiteCurrent - restrictionSitePrev;

       reference[noReferenceFragments] = new Fragment(fragmentSize, noReferenceFragments, totalSize, false);

       Fragment fragNorm = new Fragment(fragmentSize, noReferenceFragments, totalSize, false);
       Fragment fragRev = new Fragment(fragmentSize, noReferenceFragments, totalSize, true);          
       listReferenceSeeds.add(fragNorm);
       listReferenceSeeds.add(fragRev);
       noReferenceSeeds += 2;
       noReferenceFragments++;
       stepIndex++;

       totalSize =+ fragmentSize;

       // Virtual separation fragment:
       reference[noReferenceFragments] = new Fragment(RestrictionMap.ENDOFCHROMOSOMES, noReferenceFragments, -1, false);
       noReferenceFragments++;

       System.err.println("Number of fragments scaffold/chromosome " + indexChromosome++ + ": " + stepIndex);

    } // END FOR All Chromosomes


    Collections.sort(listReferenceSeeds);

    listReferenceSeedsToArray = new Fragment[listReferenceSeeds.size()];
    listReferenceSeeds.toArray(listReferenceSeedsToArray);

    TOTALFRAGMENTS = noReferenceFragments;

    System.err.println("Total number of fragments: " + TOTALFRAGMENTS);

  }

}
