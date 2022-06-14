/*
  Copyright (C) 2017 Julian Späth
  ----------------------------------------------------------------------------
  This file is part of AMMSoCK.

  AMMSoCK is free software: you can redistribute it and/or modify it under the
  terms of the GNU General Public License as published by the Free Software
  Foundation, either version 3 of the License, or (at your option) any later
  version. AMMSoCK is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
  FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
  details. You should have received a copy of the GNU General Public License
  along with AMMSoCK. If not, see <https://www.gnu.org/licenses/>.
*/

#include <algorithm>
#include <assert.h>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <time.h>
#include <vector>

using namespace std;

bool checkSpecies(string token,
                  vector<string> species) { // Prüft ob der die Spezies "Token"
                                            // in den übergebenen String vokommt
                                            // und gibt True oder False zurück
  for (int i = 0; i < species.size(); i++) {
    if (species[i].compare(token) == 0) {
      return true;
    };
  }
  return false;
}

int getIndexOfToken(
    string token,
    vector<string>
        species) { // Gibt den Index im Vektor "spezies" zurück der der Spezies
                   // "token" entspricht, falls nicht vorhanden dann -1
  for (int i = 0; i < species.size(); i++) {
    if (species[i].compare(token) == 0) {
      return i;
    };
  }
  return -1;
}

int thirdBodyOccurs(vector<int> bodyPos, int index) {
  for (int i = 0; i < bodyPos.size(); ++i) {
    if (bodyPos[i] == (index + 1)) {
      return i;
    }
  }
  return -1;
}

void getConservationMatrixAndMolarWeights(vector<string> species,
                                          vector<string> atoms,
                                          vector<double> Ma, vector<double> &Ms,
                                          vector<vector<double>> &consMatrix) {
  int ind = -1, indbeg = -1;
  // get multiplies of each atom c(a_i,s)
  for (int s = 0; s < species.size(); ++s) {
    string token = species[s];

    // delete singlet term (S)
    size_t posSinglet = token.find_first_of("()");
    if (posSinglet != -1) {
      token.erase(posSinglet);
    }

    while (token.length() > 0) {
      if (token.length() > 1) {
        // check if first alphabetical is an atom
        ind = getIndexOfToken(string(1, token[0]), atoms);
        size_t posToken = token.find_first_of("123456789");

        if (posToken == std::string::npos) {
          ind = getIndexOfToken(token, atoms);
          // if token is not found, just the first char is an atom
          if (ind == -1) {
            ind = getIndexOfToken(string(1, token[0]), atoms);
            consMatrix[ind][s] = consMatrix[ind][s] + 1.0;
            token.erase(0, 1);
          } else {
            consMatrix[ind][s] = consMatrix[ind][s] + 1.0;
            token = "";
          }
        } else {

          ind = getIndexOfToken(token.substr(0, posToken), atoms);
          // if token[0] is a species
          if (ind == -1) {
            indbeg = getIndexOfToken(string(1, token[0]), atoms);
            consMatrix[indbeg][s] = consMatrix[indbeg][s] + 1.0;
            token.erase(0, 1);
            continue;
          }
          token.erase(0, posToken);

          // check if there is there another atom
          size_t posAlpha = token.find_first_of("ABCDEFGHIJKLMNOPQRSTUVWXYZ");
          string numberOfAtom;
          if (posAlpha == std::string::npos) {
            numberOfAtom = token;
          } else {
            numberOfAtom = token.substr(0, posAlpha);
          }
          consMatrix[ind][s] =
              consMatrix[ind][s] + stof(numberOfAtom); // HIER HIER HIER//
          token.erase(0, posAlpha);
        }
      } else {
        ind = getIndexOfToken(token, atoms);
        consMatrix[ind][s] = consMatrix[ind][s] + 1.0;
        token.erase(0, 1);
      }
    }
  }

  // calculate molar weights for each species avoiding multiplication with zero
  for (int s = 0; s < species.size(); ++s) {
    for (int a = 0; a < atoms.size(); ++a) {
      if (consMatrix[a][s] > 0) {
        Ms[s] = Ms[s] + consMatrix[a][s] * Ma[a];
      }
    }
  }

  // inset mA/mS values avoiding multiplication with one
  for (int s = 0; s < species.size(); ++s) {
    for (int a = 0; a < atoms.size(); ++a) {
      consMatrix[a][s] *= Ma[a] / Ms[s];
    }
  }
}

int main(int argc, char *argv[]) {

  // check arguments
  assert(argc == 3);

  // perliminaries
  int nreac = 0, nspec = 0, natom = 0, ntb = 0;
  bool isMech = 0, isCollision = 0;
  vector<double> A, b, Ea;
  vector<string> species, atoms;
  vector<vector<double>> bodyCoeff;
  vector<int> bodyPos;
  vector<int> bodyNumber;
  vector<int> forwardReactions;
  string line, nextline, reaction, constants, reag, prod, token1, token2,
      token3, token4, token5, token6, strA, strB, strEa, colCoeffLine;

  // create streams
  ifstream mech, reacInput, molarMass;
  string path = argv[1];
  string filename = argv[2];
  mech.open((path + "/" + filename).c_str(), ios_base::in);
  ofstream reacFile("mech/data/reaction.dat");

  // get format: HOMREA *.mch or chemkin-II *.inp
  std::size_t found = filename.find_first_of(".");
  string format = filename.substr(found + 1);

  // check if mech file is available
  if (!mech) {
    cerr << "Error in parser: Source file not found.\n";
    exit(1);
  }

  if (format.compare("inp") == 0) {
    bool nextLine = true, trdBody = false, troe = false;

    while (getline(mech, line)) {

      // skip comments
      if (line[0] == '!') {
        continue;
      }

      // get elements of mechanism
      if (line.find("ELEMENTS") != std::string::npos) {
        while (getline(mech, line) && (line.find("END") == std::string::npos)) {
          stringstream str(line);
          string buffer;
          while (str >> buffer) {
            if (buffer.length() > 0) {
              atoms.push_back(buffer);
              natom++;
            }
          }
        }
        continue;
      }

      // get species of mechanism
      if (line.find("SPECIES") != std::string::npos) {
        while (getline(mech, line) && (line.find("END") == std::string::npos)) {
          stringstream str(line);
          string buffer;
          while (str >> buffer) {
            if (buffer.length() > 0) {
              species.push_back(buffer);
              nspec++;
            }
          }
        }
        continue;
      }

      // get reactions of mechanism and write them to a temporary file
      if (line.find("REACTIONS") != std::string::npos) {
        while (1) {
          // check if next line should be read
          if (nextLine) {
            if (!getline(mech, line)) {
              break;
            }
          } else {
            nextLine = true;
          }

          // break, if END is found
          if (line.find("END") != std::string::npos) {
            break;
          }

          // skip lines with "DUP" or DUPLICATE
          if ((line.find("DUP") != std::string::npos) ||
              (line.find("DUPLICATE") != std::string::npos)) {
            continue;
          }

          // check if it is an equilibrium reaction or just a forward reaction
          if ((line.find("<") == std::string::npos) &&
              (line.find("=>") != std::string::npos)) {
            forwardReactions.push_back(nreac);
            size_t posEq = line.find("=");
            line.insert(posEq, " ");
          }

          // replace '<=>' with ' = '
          std::replace(line.begin(), line.end(), '<',
                       ' '); // replace < by whitespace
          std::replace(line.begin(), line.end(), '>',
                       ' '); // replace > by whitespace

          // delete comments beginning with '!'
          size_t posComment = line.find("!");
          if (posComment != std::string::npos) {
            line.erase(posComment);
          }

          // check ocurrence of third bodies
          if (line.find("M") != std::string::npos) {
            trdBody = true;
          }

          // check if it is a troe reaction
          if (line.find("(+M)") != std::string::npos) {
            troe = true;
            // delete first brackets
            size_t posM = line.find("(+M)");
            line.erase(posM, 1);
            line.erase(posM + 2, 1);
            // delete second brackets
            posM = line.find("(+M)");
            if (posM != std::string::npos) {
              line.erase(posM, 1);
              line.erase(posM + 2, 1);
            }
          }

          // tokenize the remaining line: 'LHS' '=' 'RHS' 'A' 'b' 'Ea'
          stringstream str(line);
          str >> token1;
          str >> token2;
          str >> token3;
          str >> token4;
          str >> token5;
          str >> token6;

          // write reaction
          reacFile << token1 << " " << token2 << " " << token3 << endl;
          nreac++;
          // get A, b and Ea
          A.push_back(atof(token4.c_str()));
          b.push_back(atof(token5.c_str()));
          Ea.push_back(atof(token6.c_str()) * 4.184); // convert from cal to J

          if (troe) {
            ntb++;
            // get next line with LOW
            getline(mech, line);

            // get next line with TROE or third body coefficients
            getline(mech, line);
            size_t posTroe = line.find("TROE");
            if (posTroe != std::string::npos) {
              getline(mech, line);
            }

            // set new third body
            bodyPos.push_back(nreac);
            vector<double> efficiencies(nspec, 1.0);

            // get efficienies
            istringstream eff(line);
            string tb_species, tb_eff;
            while (getline(eff, tb_species, '/')) {
              // get next token
              if (!getline(eff, tb_eff, '/')) {
                break;
              }
              // delete whitespaces
              tb_species.erase(
                  remove(tb_species.begin(), tb_species.end(), ' '),
                  tb_species.end());
              tb_eff.erase(remove(tb_eff.begin(), tb_eff.end(), ' '),
                           tb_eff.end());

              int ind = getIndexOfToken(tb_species, species);
              if (ind != -1) {
                efficiencies[ind] = atof(tb_eff.c_str());
              } else {
                cerr << "Error in Third Body Creation: Species not found."
                     << endl;
              }
            }

            // check if current efficiencies are in bodyCoeff
            bool isContained = true;
            int ind = -1;
            for (int i = 0; i < bodyCoeff.size(); ++i) {
              for (int s = 0; s < nspec; ++s) {
                if (bodyCoeff[i][s] != efficiencies[s]) {
                  isContained = false;
                }
              }
              if (isContained) {
                ind = i;
                break;
              } else {
                isContained = true;
              }
            }
            if (ind != -1) {
              bodyNumber.push_back(ind + 1);
              ntb--;
            } else {
              bodyNumber.push_back(ntb);
              bodyCoeff.push_back(efficiencies);
            }
            nextLine = true;
            troe = false;
            trdBody = false;
            continue;
          }

          if (trdBody) {
            ntb++;
            getline(mech, nextline);

            size_t posSlash = nextline.find("/");
            if (posSlash != std::string::npos) {
              // set new third body
              bodyPos.push_back(nreac);
              vector<double> efficiencies(nspec, 1.0);

              // get efficienies
              istringstream eff(nextline);
              string tb_species, tb_eff;
              while (getline(eff, tb_species, '/')) {
                // get next token
                if (!getline(eff, tb_eff, '/')) {
                  break;
                }
                // delete whitespaces
                tb_species.erase(
                    remove(tb_species.begin(), tb_species.end(), ' '),
                    tb_species.end());
                tb_eff.erase(remove(tb_eff.begin(), tb_eff.end(), ' '),
                             tb_eff.end());

                // get index of species and store efficienies
                int ind = getIndexOfToken(tb_species, species);
                if (ind != -1) {
                  efficiencies[ind] = atof(tb_eff.c_str());
                } else {
                  cerr << "Error in Third Body Creation: Species not found."
                       << endl;
                }
              }

              // check if current efficiencies are in bodyCoeff
              bool isContained = true;
              int ind = -1;
              for (int i = 0; i < bodyCoeff.size(); ++i) {
                for (int s = 0; s < nspec; ++s) {
                  if (bodyCoeff[i][s] != efficiencies[s]) {
                    isContained = false;
                  }
                }
                // is this line equal
                if (isContained) {
                  ind = i;
                  break;
                } else {
                  isContained = true;
                }
              }
              if (ind != -1) {
                bodyNumber.push_back(ind + 1);
                ntb--;
              } else {
                bodyNumber.push_back(ntb);
                bodyCoeff.push_back(efficiencies);
              }
              nextLine = true;
            } else {

              // unspecified third body: set all efficienies to 1.0
              vector<double> efficiencies(nspec, 1.0);
              bodyPos.push_back(nreac);

              // check if current efficiencies are in bodyCoeff
              bool isContained = true;
              int ind = -1;
              for (int i = 0; i < bodyCoeff.size(); ++i) {
                for (int s = 0; s < nspec; ++s) {
                  if (bodyCoeff[i][s] != efficiencies[s]) {
                    isContained = false;
                  }
                }
                if (isContained) {
                  ind = i;
                  break;
                } else {
                  isContained = true;
                }
              }

              if (ind != -1) {
                bodyNumber.push_back(ind + 1);
                ntb--;
              } else {
                bodyNumber.push_back(ntb);
                bodyCoeff.push_back(efficiencies);
              }

              nextLine = false;
              line = nextline;
            }

            trdBody = false;
            continue;
          } // third body
        }
        continue;
      }
    } // while

  } else if (format.compare("mch") == 0) {
    // handle each line
    while (getline(mech, line)) {
      // skip duplicate reaction line
      if ((line.find("BEGIN DUPLICATE")) != -1) {
        continue;
      }

      if ((line[0] != '*') && (line[0] != '-') && (line[0] != ' ') &&
          (line[0] != 'E') && (line[0] != '0')) {
        // check if MECH block oder COLLISION block
        if (line.find("MECHANISM") != -1) {
          isMech = 1;
          continue;
        }
        if (line.find("COLLISION") != -1) {
          isMech = 0;
          isCollision = 1;
          continue;
        }

        // handle MECHANISM block
        if (isMech) {
          nreac++;
          reaction = line.substr(0, 44);

          // check if it is only a forward reaction
          if (reaction.find(">") != std::string::npos) {
            forwardReactions.push_back(nreac - 1);
          }
          std::replace(reaction.begin(), reaction.end(), '>',
                       '='); // replace > by =

          // handle reaction
          token1 = reaction.substr(0, 8);
          token2 = reaction.substr(9, 8);
          token3 = reaction.substr(18, 8);
          token4 = reaction.substr(27, 8);
          token5 = reaction.substr(36, 8);

          // delete white spaces
          token1.erase(remove(token1.begin(), token1.end(), ' '), token1.end());
          token2.erase(remove(token2.begin(), token2.end(), ' '), token2.end());
          token3.erase(remove(token3.begin(), token3.end(), ' '), token3.end());
          token4.erase(remove(token4.begin(), token4.end(), ' '), token4.end());
          token5.erase(remove(token5.begin(), token5.end(), ' '), token5.end());

          if (!checkSpecies(token1, species) && (token1[0] != 'M') &&
              (token1.size() > 0)) {
            species.push_back(token1);
          }
          if (!checkSpecies(token2, species) && (token2[0] != 'M') &&
              (token2.size() > 0)) {
            species.push_back(token2);
          }
          if (!checkSpecies(token3, species) && (token3[0] != 'M') &&
              (token3.size() > 0)) {
            species.push_back(token3);
          }
          if (!checkSpecies(token4, species) && (token4[0] != 'M') &&
              (token4.size() > 0)) {
            species.push_back(token4);
          }
          if (!checkSpecies(token5, species) && (token5[0] != 'M') &&
              (token5.size() > 0)) {
            species.push_back(token5);
          }

          // handle reaction constants
          constants = line.substr(46);
          strA = constants.substr(0, 9);
          strB = constants.substr(10, 17);
          strEa = constants.substr(17);
          A.push_back(atof(strA.c_str()));
          b.push_back(atof(strB.c_str()));
          Ea.push_back(atof(strEa.c_str()) * 1e+3); // convert from kJ to J

          // look for third bodys and store number as well as position
          size_t posThird = reaction.find("M");
          if (posThird != std::string::npos) {
            bodyPos.push_back(nreac);
            bodyNumber.push_back((int)reaction[posThird + 2] - '0');
            reaction.erase(posThird + 1, 3);
            posThird = reaction.find("M", posThird + 2);
            if (posThird != std::string::npos) {
              reaction.erase(posThird + 1, 3);
            }
          }

          stringstream str(reaction);
          string newreac = "";
          while (str >> token1) {
            newreac += token1;
          }

          size_t posEq = newreac.find("=");
          if (posEq != std::string::npos) {
            newreac.insert(posEq, " ");
            newreac.insert(posEq + 2, " ");
          }

          reacFile << newreac << endl;

        } // isMech

        // get number of species within reaction
        nspec = species.size();

        // handle COLLISION block
        if (isCollision) {
          // get collision efficiencies and their quantity
          int ncolCoeff = (line.size() - 9) / 9 + 1;
          getline(mech, colCoeffLine);

          vector<double> vTemp(nspec, 1.0);
          map<string, double> colCoeff;
          map<string, double>::iterator it;

          for (int i = 1; i <= ncolCoeff; ++i) {
            // get tokens
            token1 = line.substr(9 * i, 8);
            token2 = colCoeffLine.substr(9 * i, 8);

            // delete whitespaces
            token1.erase(remove(token1.begin(), token1.end(), ' '),
                         token1.end());
            token2.erase(remove(token2.begin(), token2.end(), ' '),
                         token2.end());

            // create new map and insert into bodyCoeff
            for (int k = 0; k < nspec; ++k) {
              if ((token1.compare(species[k])) == 0) {
                vTemp[k] = atof(token2.c_str());
              }
            }
          }

          bodyCoeff.push_back(vTemp);
        } // isCollision
      }
    } // while

    // get atoms: WARNING! Atoms with more than two chars, e.g. 'He' or 'Ar'
    // have to have non-capital chars. HE is not valid for helium.
    for (int s = 0; s < nspec; ++s) {
      token1 = species[s];
      while (token1.length() > 0) {
        if ((token1[0] != '0') && (token1[0] != '1') && (token1[0] != '2') &&
            (token1[0] != '3') && (token1[0] != '4') && (token1[0] != '5') &&
            (token1[0] != '6') && (token1[0] != '7') && (token1[0] != '8') &&
            (token1[0] != '9')) {
          size_t posNonCapital =
              token1.find_first_of("abcdefghijklmnopqrstuvwxyz");
          string atom;
          if (posNonCapital != std::string::npos) {
            atom = token1.substr(0, posNonCapital + 1);
            token1.erase(0, posNonCapital + 1);
          } else {
            atom = token1.substr(0, 1);
            token1.erase(0, 1);
          }
          if (!checkSpecies(atom, atoms)) {
            atoms.push_back(atom);
          }
        } else {
          token1.erase(0, 1);
        }
      }
    }

  } else {
    cerr << "Error in parser: Unkown format. Please use HOMREA (*.mch) or "
            "chemkin-II (*.inp) format.\n";
    exit(1);
  }

  // assign nreac and natom
  natom = atoms.size();
  nspec = species.size();

  // begin post-processing and parse reaction file
  reacFile.close();
  nspec = species.size();
  reacInput.open("mech/data/reaction.dat", ios_base::in);
  vector<vector<int>> nuprime;
  vector<vector<int>> nu2prime;
  for (int i = 0; i < nspec; ++i) {
    nuprime.push_back(vector<int>(nreac, 0));
    nu2prime.push_back(vector<int>(nreac, 0));
  }

  for (int r = 1; r <= nreac; ++r) {
    // tokenize line
    getline(reacInput, line);
    stringstream str(line);
    str >> token1; // left-hand side
    str >> token2; // equality
    str >> token3; // right-hand side

    // tokenize and deal with left-hand side
    vector<string> tokens;
    string curToken;
    istringstream istr_lhs(token1);
    int prefactor = 0, molecularity = 0, ind = -1;

    while (getline(istr_lhs, token4, '+')) {
      tokens.push_back(token4);
    }

    for (int i = 0; i < tokens.size(); ++i) {
      prefactor = 0;
      curToken = tokens[i];

      // check first char - note: higher than molecularity 3 is nearly
      // impossible
      if ((curToken[0] == '2') || (curToken[0] == '3') ||
          (curToken[0] == '4') || (curToken[0] == '5') ||
          (curToken[0] == '6') || (curToken[0] == '7') ||
          (curToken[0] == '8') || (curToken[0] == '9')) {
        prefactor = (int)curToken[0] - '0';
        curToken.erase(0, 1);
      } else {
        prefactor = 1;
      }

      // find index of species
      ind = getIndexOfToken(curToken, species);
      // ignore third bodies and set nuprime
      if (ind != -1) {
        nuprime[ind][r - 1] += prefactor;
      }

      // assign molecularity
      molecularity += prefactor;
    }

    // tokenize and deal with right-hand side
    istringstream istr_rhs(token3);
    tokens.clear();

    while (getline(istr_rhs, token4, '+')) {
      tokens.push_back(token4);
    }

    for (int i = 0; i < tokens.size(); ++i) {
      prefactor = 0;
      curToken = tokens[i];

      // check first char - note: higher than molecularity 3 is nearly
      // impossible
      if ((curToken[0] == '2') || (curToken[0] == '3') ||
          (curToken[0] == '4') || (curToken[0] == '5') ||
          (curToken[0] == '6') || (curToken[0] == '7') ||
          (curToken[0] == '8') || (curToken[0] == '9')) {
        prefactor = (int)curToken[0] - '0';
        curToken.erase(0, 1);
      } else {
        prefactor = 1;
      }

      // find index of species
      ind = getIndexOfToken(curToken, species);
      // ignore third bodies and set nu2prime
      if (ind != -1) {
        nu2prime[ind][r - 1] += prefactor;
      }
    }

    // adjust preexponential factor due to units: convert from cm^3 to m^3
    if (molecularity == 2) {
      A[r - 1] *= 1e-6;
    } else if (molecularity == 3) {
      A[r - 1] *= 1e-12;
    }
  }
  reacInput.close();

  // get molar masses
  vector<double> Ms(nspec, 0.0);
  vector<double> Ma(natom, 0.0);
  int ind = -1;
  line = "";

  molarMass.open((path + "/molarWeights.dat").c_str(), ios_base::in);
  while (getline(molarMass, line)) {
    // skip comment lines
    if (line[0] == '#') {
      continue;
    }

    // handle line
    istringstream iss(line);
    iss >> token1;
    iss >> token2;

    ind = getIndexOfToken(token1, atoms);
    if (ind != -1) {
      Ma[ind] = atof(token2.c_str());
    }
  }

  // get Conservation Matrix
  vector<vector<double>> consMatrix(natom, vector<double>(nspec, 0.0));
  getConservationMatrixAndMolarWeights(species, atoms, Ma, Ms, consMatrix);

  // open thermo.dat and get NASA coefficients
  ifstream thermo;
  thermo.open((path + "/thermo.dat").c_str(), ios_base::in);

  vector<vector<double>> nasa(nspec, vector<double>(14, 0.0));
  vector<double> switchingPoints(nspec, 0.0);

  while (getline(thermo, line)) {

    // skip other lines
    if (line.length() != 80) {
      continue;
    }

    // check first line and get species and switching points
    if (line[79] == '1') {
      // get specie and delete whitespaces
      int specInd = -1;
      token1 = line.substr(0, 15);
      token1.erase(remove(token1.begin(), token1.end(), ' '), token1.end());
      for (int i = 0; i < nspec; ++i) {
        if (species[i].compare(token1) == 0) {
          specInd = i;
        }
      }

      // if species is not found, skip the next three lines
      if (specInd == -1) {
        getline(thermo, line);
        getline(thermo, line);
        getline(thermo, line);
        continue;
      }

      // get switching point and delete whitespaces
      token2 = line.substr(64, 10);
      token2.erase(remove(token2.begin(), token2.end(), ' '), token2.end());
      switchingPoints[specInd] = atof(token2.c_str());

      // deal with 2nd line
      getline(thermo, line);
      token1 = line.substr(0, 15);
      nasa[specInd][7] = atof(token1.c_str());
      token2 = line.substr(15, 15);
      nasa[specInd][8] = atof(token2.c_str());
      token3 = line.substr(30, 15);
      nasa[specInd][9] = atof(token3.c_str());
      token4 = line.substr(45, 15);
      nasa[specInd][10] = atof(token4.c_str());
      token5 = line.substr(60, 15);
      nasa[specInd][11] = atof(token5.c_str());

      // deal with 3rd line
      getline(thermo, line);
      token1 = line.substr(0, 15);
      nasa[specInd][12] = atof(token1.c_str());
      token2 = line.substr(15, 15);
      nasa[specInd][13] = atof(token2.c_str());
      token3 = line.substr(30, 15);
      nasa[specInd][0] = atof(token3.c_str());
      token4 = line.substr(45, 15);
      nasa[specInd][1] = atof(token4.c_str());
      token5 = line.substr(60, 15);
      nasa[specInd][2] = atof(token5.c_str());

      // deal with 4th line
      getline(thermo, line);
      token1 = line.substr(0, 15);
      nasa[specInd][3] = atof(token1.c_str());
      token2 = line.substr(15, 15);
      nasa[specInd][4] = atof(token2.c_str());
      token3 = line.substr(30, 15);
      nasa[specInd][5] = atof(token3.c_str());
      token4 = line.substr(45, 15);
      nasa[specInd][6] = atof(token4.c_str());
    }
  } // end get lines from thermo
  thermo.close();

  // write mech data to files
  ofstream mechdata("mech/data/mechdata.dat");

  mechdata << nspec << endl;
  for (int i = 0; i < nspec; ++i) {
    mechdata << setprecision(15) << species[i] << " " << Ms[i] << endl;
  }

  mechdata << natom << endl;
  for (int i = 0; i < natom; ++i) {
    mechdata << atoms[i] << " " << Ma[i] << endl;
  }

  mechdata << nreac << endl;
  for (int i = 0; i < nuprime.size(); ++i) {
    for (int k = 0; k < nuprime[0].size(); ++k) {
      mechdata << nuprime[i][k] << " ";
    }
    mechdata << endl;
  }
  for (int i = 0; i < nu2prime.size(); ++i) {
    for (int k = 0; k < nu2prime[0].size(); ++k) {
      mechdata << nu2prime[i][k] << " ";
    }
    mechdata << endl;
  }

  mechdata << ntb << endl;
  for (int i = 0; i < bodyCoeff.size(); ++i) {
    for (int k = 0; k < bodyCoeff[0].size(); ++k) {
      mechdata << bodyCoeff[i][k] << " ";
    }
    mechdata << endl;
  }

  mechdata << bodyPos.size() << endl;
  for (int i = 0; i < bodyPos.size(); ++i) {
    mechdata << bodyPos[i] << " " << bodyNumber[i] << endl;
  }

  mechdata << forwardReactions.size() << endl;
  for (int i = 0; i < forwardReactions.size(); ++i) {
    mechdata << setprecision(15) << forwardReactions[i] << endl;
  }

  for (int i = 0; i < natom; ++i) {
    for (int k = 0; k < nspec; ++k) {
      mechdata << setprecision(15) << consMatrix[i][k] << " ";
    }
    mechdata << endl;
  }
  mechdata.close();

  // write thermo data to files
  ofstream thermodata("mech/data/thermodata.dat");

  for (int i = 0; i < A.size(); ++i) {
    thermodata << std::setprecision(15) << A[i] << " " << b[i] << " " << Ea[i]
               << endl;
  }
  for (int i = 0; i < nspec; ++i) {
    for (int k = 0; k < 14; ++k) {
      thermodata << nasa[i][k] << " ";
    }
    thermodata << switchingPoints[i] << endl;
  }
  thermodata.close();

  cout << "\n*** readData: Work done. ***" << endl << endl;
  cout << "Statistics " << endl;
  cout << "\t Mechanism: " << filename << endl;
  cout << "\t Reactions: " << nreac << endl;
  cout << "\t #Species: " << nspec << endl;
  cout << "\t Species: ";
  for (int i = 0; i < nspec; ++i) {
    cout << species[i] << " ";
  }
  cout << endl;
  cout << "\t #Atoms: " << natom << endl;
  cout << "\t Atoms: ";
  for (int i = 0; i < natom; ++i) {
    cout << atoms[i] << " ";
  }
  cout << endl;
  cout << "\t #Third Bodies: " << bodyCoeff.size() << endl;

  mech.close();

  // write reader statistics to file
  ofstream stats("mech/data/readerStats.dat");
  stats << nreac << endl;
  stats << nspec << endl;
  for (int i = 0; i < nspec; ++i) {
    stats << species[i] << endl;
  }
  stats << natom << endl;
  for (int i = 0; i < natom; ++i) {
    stats << atoms[i] << endl;
  }
  stats << ntb << endl;
  stats << bodyPos.size() << endl;
  stats << forwardReactions.size() << endl;
  for (int i = 0; i < forwardReactions.size(); ++i) {
    stats << forwardReactions[i] << endl;
  }
  stats.close();

  return 0;
}
