#include "zone.h"␊
␊
#include "district.h"␊
#include "building.h"␊
#include "models.h"␊
#include <algorithm>␊
␊
#include <array>␊
#include <cmath>␊
#include <limits>␊
#include "occupants.h"␊
#include "scene.h"␊
␊
// *** Zone class, CitySim   *** //␊
// *** jerome.kaempf@epfl.ch *** //␊
␊
Zone::Zone(unsigned int id, Building* pBuilding, bool groundFloor, float Vi, vector<Wall*> walls, vector<Roof*> roofs, vector<Surface*> surfaces, vector<Floor*> floors, Occupants *pOccupants)␊
    :id(id),pBuilding(pBuilding),walls(walls),roofs(roofs),surfaces(surfaces),floors(floors),Vi(Vi),groundFloor(groundFloor),pOccupants(pOccupants),logStream(pBuilding->logStream.rdbuf()) {␊
␊
    logStream << "Zone creator." << endl << flush;␊
␊
    // the zone cannot have a volume of zero␊
    if (Vi<=0.f) throw(string("In Zone id ") + toString(id) + string(": Cannot create a zone without a positive volume value."));␊
␊
␊
    for (unsigned int i=0; i<walls.size(); ++i) {␊
        walls[i]->setBuildingRef(pBuilding);␊
    }␊
    for (unsigned int i=0; i<roofs.size(); ++i) {␊
        roofs[i]->setBuildingRef(pBuilding);␊
    }␊
    for (unsigned int i=0; i<floors.size(); ++i) {␊
        floors[i]->setBuildingRef(pBuilding);␊
    }␊
    for (unsigned int i=0; i<surfaces.size(); ++i) {␊
        surfaces[i]->setBuildingRef(pBuilding);␊
    }␊
␊
    update(true);␊
␊
}␊
␊
/**␊
 * @brief Zone::updateSimulationModelParameters: update all "secondary" parameters, i.e.␊
 *  parameters which are based on primary parameters that can be modified in the interface CitySimPro␊
 */␊
void Zone::update(bool constructor){␊
    // note: the computation of Qsun1 and Qsun2 is made when a call to the method getQsun1 and getQsun2 is made␊
    if(!constructor){␊
        Swa=0.;␊
        Swi=0.;␊
        SwiO=0.;␊
        Sro=0.f;␊
        Kroof = 0.f;␊
        Kwindow=0.f;␊
    }␊
    // loop on the surfaces to get the surface of walls (Swa) and the surface of windows (Swi)␊
    for (unsigned int i=0; i<walls.size(); ++i) {␊
        logStream << "Wall #" << i << " Glazing Ratio: " << walls[i]->getGlazingRatio() << "\tArea: " << walls[i]->getArea() << endl << flush;␊
        Swa  += walls[i]->getWallArea();␊
        Swi  += walls[i]->getGlazingArea();␊
        SwiO += walls[i]->getGlazingArea()*walls[i]->getGlazingOpenableRatio();␊
        Kwindow += walls[i]->getGlazingUvalue()*walls[i]->getGlazingArea();␊
        walls[i]->updateGlazingGvalueHemispherical();␊
    }␊
    for (unsigned int i=0; i<roofs.size(); ++i) {␊
        logStream << "Roof #" << i << " Glazing Ratio: " << roofs[i]->getGlazingRatio() << "\tArea: " << roofs[i]->getArea() << endl << flush;␊
        Sro += roofs[i]->getRoofArea();␊
        Swi += roofs[i]->getGlazingArea();␊
        SwiO += roofs[i]->getGlazingArea()*roofs[i]->getGlazingOpenableRatio();␊
        Kroof += roofs[i]->getComposite()->getUvalue()*roofs[i]->getRoofArea();␊
        Kwindow += roofs[i]->getGlazingUvalue()*roofs[i]->getGlazingArea();␊
        roofs[i]->updateGlazingGvalueHemispherical();␊
    }␊
    logStream << "Total wall surface (m^2): " << Swa << "\tTotal window surface (m^2): " << Swi << "\tTotal roof surface (m^2): " << Sro << endl << flush;␊
    logStream << "Total window conductance (W/K): " << Kwindow << "\tTotal roof conductance (W/K): " << Kroof << endl << flush;␊
␊
    // computation of Ww and Wa, as a first approximation, the ratio of the surfaces to the total surface, JK - 18.04.2009␊
    //if ( (Swa+Swi) > 0. ) { Ww=Swa/(Swa+Swi); Wa=Swi/(Swa+Swi); }␊
    //logStream << "Ww: " << Ww << "\tWa: " << Wa << endl << flush;␊
␊
    if (pBuilding->getDistrict()->getScene()->getClimate()!=NULL)␊
        setAirDensity(pBuilding->getDistrict()->getScene()->getClimate()->getAirDensity());␊
    else␊
        setAirDensity(Climate::getAirDensity(0.f)); // if no Climate make an assumption that the altitude of the situation is 0 m (sea level)␊
␊
    // output of the values␊
    logStream << "Vi: " << Vi << endl << flush;␊
}␊
␊
vector<Surface*> Zone::getAllSurfaces(){␊
    vector<Surface*> v = surfaces;␊
    v.insert(v.begin(),walls.begin(),walls.end());␊
    v.insert(v.begin(),floors.begin(),floors.end());␊
    v.insert(v.begin(),roofs.begin(),roofs.end());␊
    return v;␊
}␊
␊
void Zone::computeVolume() {␊
␊
    // initialize the air volume␊
    Vi = 0.f;␊
    // computes the volume under the surfaces␊
    for (vector<Wall*>::iterator it=walls.begin();it!=walls.end();++it) Vi += (*it)->getVolume();␊
    for (vector<Roof*>::iterator it=roofs.begin();it!=roofs.end();++it) Vi += (*it)->getVolume();␊
    for (vector<Floor*>::iterator it=floors.begin();it!=floors.end();++it) Vi += (*it)->getVolume();␊
␊
}␊
␊
float Zone::getTmin(unsigned int day, unsigned int hour) {␊
    try {␊
        // if a TemperatureProfile exists then return the temperature in it otherwise Tmin␊
        if (Tprofile) return pBuilding->getDistrict()->getTemperatureProfiles()->getYearProfile(*Tprofile)->getDayProfile(day)->getHourValue_Tmin(hour);␊
        else return Tmin;␊
    }␊
    catch(exception& e) {␊
        return Tmin;␊
    }␊
}␊
␊
float Zone::getTmax(unsigned int day, unsigned int hour) {␊
    try {␊
        // if a TemperatureProfile exists then return the temperature in it otherwise Tmin␊
        if (Tprofile) return pBuilding->getDistrict()->getTemperatureProfiles()->getYearProfile(*Tprofile)->getDayProfile(day)->getHourValue_Tmax(hour);␊
        else return Tmax;␊
    }␊
    catch(exception& e) {␊
        return Tmax;␊
    }␊
}␊
␊
double Zone::getKappa3() {␊
␊
    double kappa3 = 0.;␊
    for (size_t i=0; i<roofs.size(); ++i){␊
        kappa3 += roofs[i]->getKappa()*roofs[i]->get_hr()*roofs[i]->getRoofArea();␊
        //cout << "Zone::getKappa3: Kappa="<< roofs[i]->getKappa()<<" hr="<<roofs[i]->get_hr() << " area="<<roofs[i]->getRoofArea()<<endl;␊
    }␊
    return kappa3;␊
␊
}␊
␊
␊
float Zone::getQsun1() {␊
    // computes Qsun1 (what goes on the enveloppe of the building) from the values on the facades␊
    float value1 = 0.f;␊
    for (size_t i=0; i<walls.size(); ++i) { // loop on the external walls␊
        value1 += ( walls[i]->getWallArea()␊
                    *(walls[i]->getShortWaveIrradiance()*(1.-walls[i]->getShortWaveReflectance_opaque())␊
                      +walls[i]->get_hr()*walls[i]->getEnvironmentalTemperature()) );␊
    }␊
    return value1;␊
}␊
␊
float Zone::getQsun2() {␊
    // computes Qsun2 (what goes inside the building with shading devices in use) from the values on the facades␊
    float value2 = 0.f;␊
    for (size_t i=0; i<walls.size(); ++i) { // loop on the external walls␊
        //logStream << "Irradiation: " << walls[i]->getShortWaveIrradiation()/walls[i]->getArea()␊
        value2 += ( walls[i]->getGlazingRatio()*walls[i]->getArea()*␊
                    (walls[i]->getGlazingGvalue(walls[i]->getBeamAngle())*walls[i]->getBeamIrradiance()␊
                     +walls[i]->getGlazingGvalueHemispherical()*(walls[i]->getShortWaveIrradiance()-walls[i]->getBeamIrradiance()))*␊
                    (0.15f+0.85f*walls[i]->getLowerShadingState()) );␊
␊
        //cout << "Wall: " << walls[i]->getId() << " Irradiance: " << walls[i]->getShortWaveIrradiation()/walls[i]->getArea() << " shading state: " << walls[i]->getLowerShadingState() << "\tshading: " << (0.15f+0.85f*walls[i]->getLowerShadingState()) << endl;␊
    }␊
    for (size_t i=0; i<roofs.size(); ++i) { // loop on the external roofs␊
        value2 += ( roofs[i]->getGlazingRatio()*roofs[i]->getArea()*␊
                    (roofs[i]->getGlazingGvalue(roofs[i]->getBeamAngle())*roofs[i]->getBeamIrradiance()␊
                     +roofs[i]->getGlazingGvalueHemispherical()*(roofs[i]->getShortWaveIrradiance()-roofs[i]->getBeamIrradiance()))*␊
                    (0.15f+0.85f*roofs[i]->getShadingState()) );␊
␊
        //cout << "Roof: " << roofs[i]->getId() << " Irradiance: " << roofs[i]->getShortWaveIrradiation()/roofs[i]->getArea() << " shading state: " << roofs[i]->getShadingState() << "\tshading: " << (0.15f+0.85f*roofs[i]->getShadingState()) << endl;␊
    }␊
    return value2;␊
}␊
␊
float Zone::getQsun3() {␊
    // computes Qsun3 (what goes on the envelope of the building) from the values on the roofs␊
    float Qsun3 = 0.f;␊
    for (size_t i=0; i<roofs.size(); ++i) {␊
        Qsun3 += roofs[i]->getKappa()*roofs[i]->getRoofArea()␊
                 *( roofs[i]->getShortWaveIrradiance()*(1.-roofs[i]->getShortWaveReflectance_opaque())␊
                    +roofs[i]->get_hr()*roofs[i]->getEnvironmentalTemperature() );␊
    }␊
    return Qsun3;␊
}␊
␊
double Zone::getTa(unsigned int day, unsigned int hour) {␊
    // pre-simulation results are not erased, but part of the vector might have been saved already␊
    return Ta.at((day-1)*24 + hour -1 + pBuilding->getDistrict()->getScene()->getPreTimeStepsSimulated() - pBuilding->getDistrict()->getScene()->getSimulationIndex());␊
}␊
␊
void Zone::setOccupantsCountAndActivity(unsigned int day, unsigned int hour)␊
{␊
    // initialises the sensible and radiative gains Lc and Lr␊
    Lc = 0.f;␊
    Lr = 0.f;␊
␊
    // two different loops, if stochastic or not␊
    if (occupantsStochastic) {␊
        // loop to determine the occupants' count␊
        occupantsCount = 0.f; // has to be a float due to deterministic procedures␊
        for (unsigned int i=0; i<static_cast<unsigned int>(round(occupantsNumber)); ++i) {␊
            if (randomUniform(0,1) <= occupantsYearProfile->getDayProfile(day)->getHourValue(hour)) occupantsCount+=1.f;␊
        }␊
        // loop to determine the occupants' behavior␊
        float activity, randomNumber;␊
        unsigned int deviceType;␊
        for (size_t n=0; n < occupantsCount; ++n) {␊
            activity = 0.f;␊
            randomNumber = randomUniform(0,1); // draw one random number per occupant to determine the activity performed␊
            if (activityType != numeric_limits<unsigned int>::signaling_NaN()) { // activity is defined␊
                for(size_t j=0; j < pBuilding->getDistrict()->getActivityType(activityType)->getnActivities(); ++j) {␊
                    activity += pBuilding->getDistrict()->getActivityType(activityType)->getActivityProbability(j,hour);␊
                    if (randomNumber <= activity) {␊
                        // execute the activity␊
                        deviceType = pBuilding->getDistrict()->getActivityType(activityType)->getActivityDeviceType(j);␊
                        // loop on all devices listed in this activity␊
                        for (size_t k = 0; k < pBuilding->getDistrict()->getDeviceType(deviceType)->getnDevices(); ++k) {␊
                            if (randomUniform(0,1) <= pBuilding->getDistrict()->getDeviceType(deviceType)->getDeviceProbability(k,hour)) {␊
                                // this device is used, add its electricity consumption to the building␊
                                pBuilding->addElectricConsumption(float(Model::dt)*pBuilding->getDistrict()->getDeviceType(deviceType)->getDeviceAvgPower(k));␊
                                // adds the sensible and convective gains due to the device use␊
                                Lc += pBuilding->getDistrict()->getDeviceType(deviceType)->getDeviceAvgPower(k)*pBuilding->getDistrict()->getDeviceType(deviceType)->getDeviceConvectiveFraction(k);␊
                                Lr += pBuilding->getDistrict()->getDeviceType(deviceType)->getDeviceAvgPower(k)*pBuilding->getDistrict()->getDeviceType(deviceType)->getDeviceRadiativeFraction(k);␊
                            }␊
                        }␊
                        // once the activity is accomplished, break the loop on activities␊
                        break;␊
                    }␊
                }␊
            }␊
        }␊
␊
    }␊
    else {␊
        // define the occupants count␊
        occupantsCount = occupantsNumber * occupantsYearProfile->getDayProfile(day)->getHourValue(hour);␊
        // loop to determine occupants behaviour␊
        unsigned int deviceType;␊
        float electricity;␊
        if (activityType != numeric_limits<unsigned int>::signaling_NaN()) { // activity is defined␊
            for(size_t j=0; j < pBuilding->getDistrict()->getActivityType(activityType)->getnActivities(); ++j) {␊
                // loop on all devices listed in this activity␊
                deviceType = pBuilding->getDistrict()->getActivityType(activityType)->getActivityDeviceType(j);␊
                for (size_t k = 0; k < pBuilding->getDistrict()->getDeviceType(deviceType)->getnDevices(); ++k) {␊
                    // this device is used, add its electricity consumption to the building␊
                    electricity = occupantsCount␊
                                  *pBuilding->getDistrict()->getActivityType(activityType)->getActivityProbability(j,hour)␊
                                  *pBuilding->getDistrict()->getDeviceType(deviceType)->getDeviceProbability(k,hour)␊
                                  *pBuilding->getDistrict()->getDeviceType(deviceType)->getDeviceAvgPower(k);␊
                    pBuilding->addElectricConsumption(float(Model::dt)*electricity);␊
                    // adds the sensible and convective gains due to the device use␊
                    Lc += electricity*pBuilding->getDistrict()->getDeviceType(deviceType)->getDeviceConvectiveFraction(k);␊
                    Lr += electricity*pBuilding->getDistrict()->getDeviceType(deviceType)->getDeviceRadiativeFraction(k);␊
                }␊
            }␊
        }␊
    }␊
␊
    // sensible and radiative gains due to the occupants␊
    Lc += occupantsCount*occupantsSensibleHeat*(1.f-occupantsSensibleHeatRadiantFraction);␊
    Lr += occupantsCount*occupantsSensibleHeat*occupantsSensibleHeatRadiantFraction;␊
␊
}␊
␊
float Zone::getDHWConsumption(unsigned int day, unsigned int hour)␊
{␊
    // do some tests␊
    if (!dhwYearProfile) {␊
        throw(string("DHWYearProfile not defined for Zone id=")+toString(id)+" in Building id="+toString(this->getpBuilding()->getId())+".");␊
    }␊
␊
    if (occupantsStochastic) {␊
        float DHWconsumption = 0.f;␊
        for (unsigned int i=0; i<occupantsCount; ++i) {␊
            if (randomUniform(0,1) <= dhwYearProfile->getDayProfile(day)->getHourValue(hour)) DHWconsumption+=dhwYearProfile->getDayProfile(day)->getWaterConsumption()*dhwYearProfile->getDayProfile(day)->getHourValue(hour); // Cognet: Deleted the division by 24 and times "getHourValue", so that the expected consumed volume per person after one day is "getWaterConsumption". Basically each person probabilistically consumes the whole "getWaterConsumption" liters in one go. TODO check this.␊
        }␊
        return DHWconsumption;␊
    }␊
    else return occupantsNumber*occupantsYearProfile->getDayProfile(day)->getHourValue(hour)*dhwYearProfile->getDayProfile(day)->getWaterConsumption()*dhwYearProfile->getDayProfile(day)->getHourValue(hour); // Cognet: Deleted the division by 24, so that the probabilities are normalized. TODO: check that this is correct.␊
}␊
␊
void Zone::writeXML(ofstream& file, string tab=""){␊
␊
    file << tab << "<Zone id=\"" << id << "\" volume=\"" << Vi << "\" psi=\"" << Kpsi␊
         << "\" Tmin=\"" << Tmin << "\" Tmax=\"" << Tmax << "\" groundFloor=\"";␊
␊
    if (groundFloor)␊
        file << "true";␊
    else␊
        file << "false";␊
    file << "\" nightVentilationBegin=\"" << nightVentilationBegin << "\" nightVentilationEnd=\"" << nightVentilationEnd;␊
    file << "\">" << endl;␊
    string subtab =tab+"\t";␊
    file << subtab << "<Occupants n=\"" << occupantsNumber << "\" sensibleHeat=\"" << getOccupantsSensibleHeat()␊
                                                           << "\" sensibleHeatRadiantFraction=\"" << getOccupantsSensibleHeatRadiantFraction()␊
                                                           << "\" latentHeat=\"" << getOccupantsLatentHeat();␊
    file << "\" type=\"" << occupantsYearProfile->getId() << "\"";␊
    if (activityType != numeric_limits<unsigned int>::signaling_NaN()) file << " activityType=\"" << activityType << "\"";␊
    file << " DHWType=\"" << getDHWYearProfile()->getId() << "\"";␊
    file << "/>" << endl;␊
␊
    for (unsigned int i=0; i<walls.size(); ++i){␊
        walls[i]->writeXML(file,subtab);␊
    }␊
    for (unsigned int i=0; i<roofs.size(); ++i){␊
        roofs[i]->writeXML(file,subtab);␊
    }␊
    for (unsigned int i=0; i<floors.size(); ++i){␊
        floors[i]->writeXML(file,subtab);␊
    }␊
    // Write Shading Surfaces inside the building␊
    for (size_t i=0; i<surfaces.size(); ++i) {␊
        file << subtab << "<Surface id=\"" << surfaces[i]->getId();␊
        file << "\" ShortWaveReflectance=\"" << surfaces[i]->getShortWaveReflectance();␊
        file << "\">" << endl;␊
        surfaces[i]->writeXML(file,subtab+"\t");␊
        file << subtab << "</Surface>" << endl;␊
    }␊
␊
    file << tab << "</Zone>" << endl;␊
}␊
␊
// Zone2N␊
Zone2N::Zone2N(unsigned int id, Building* pBuilding, bool groundFloor, float Vi, vector<Wall*> walls, vector<Roof*> roofs, vector<Surface*> surfaces, vector<Floor*> floors, Occupants *pOccupants)␊
        :Zone(id,pBuilding,groundFloor,Vi,walls,roofs,surfaces,floors,pOccupants) {␊
␊
    // intialisation␊
    nNodes=2;␊
    Tw=15.f;␊
␊
    update(true);␊
␊
}␊
␊
/**␊
 * @brief Zone2N::updateSimulationModelParameters: update all "secondary" parameters, i.e.␊
 *  parameters which are based on other parameters that can be modified in the interface CitySimPro␊
 */␊
void Zone2N::update(bool constructor){␊
    if(!constructor){␊
        Zone::update();␊
    }␊
    // If constructor, this was already called when the parent Zone was initialised␊
␊
    if (walls.empty()) { // no external wall in this zone, temperature decoupled from the model␊
        Cw = 1.f;␊
        Kw1 = 1.f;␊
        Kw2 = 1.f;␊
    }␊
    else {␊
        walls[0]->getComposite()->getSimplifiedNode(Cw, Kw1, Kw2);␊
        // multiply with the total walls area␊
        Cw *= Swa;␊
    }␊
␊
    // outputs the values␊
    //logStream << "Cw: " << Cw << "\tKw1: " << Kw1*Swa << "\tKw2: " << Kw2*Swa << endl;␊
}␊
␊
double Zone2N::getKappa1() {␊
␊
    // if no walls, returns a dummy value to couple with the outdoor environment␊
    if (walls.empty()) return 1.;␊
␊
    double kappa1 = 0.;␊
    for (size_t i=0; i<walls.size(); ++i)␊
        kappa1 += Kw1*walls[i]->getWallArea()*(walls[i]->get_hc()+walls[i]->get_hr())/(Kw1+walls[i]->get_hc()+walls[i]->get_hr());␊
    return kappa1;␊
␊
}␊
␊
void Zone2N::setTos(float Tout) {␊
    // Ta.back() stores the air node temperature from the most recent thermal solve.␊
    // setTos runs immediately after the solver so roofs are coupled to the current-step␊
    // indoor air temperature rather than a lagged value.␊
    // define the surface temperature per surface␊
    for (size_t i=0; i<walls.size(); ++i) {␊
        float wallTemperature = ( Kw1*Tw + walls[i]->get_hc()*Tout␊
                                  +walls[i]->getShortWaveIrradiance()*(1.-walls[i]->getShortWaveReflectance_opaque())␊
                                  +walls[i]->get_hr()*walls[i]->getEnvironmentalTemperature() )␊
                                / (Kw1 + walls[i]->get_hc() + walls[i]->get_hr());␊
        walls[i]->setTemperature(wallTemperature);␊
    }␊
    for (size_t i=0; i<roofs.size(); ++i) {␊
        float roofTemperature = ( roofs[i]->getKr()*Ta.back() + roofs[i]->get_hc()*Tout␊
                                  +roofs[i]->getShortWaveIrradiance()*(1.-roofs[i]->getShortWaveReflectance_opaque())␊
                                  +roofs[i]->get_hr()*roofs[i]->getEnvironmentalTemperature() )␊
                                / (roofs[i]->getKr() + roofs[i]->get_hc() + roofs[i]->get_hr());␊
        roofs[i]->setTemperature(roofTemperature);␊
    }␊
}␊
␊
float Zone2N::computeBipvHeatingGain(Wall& wall, float Tout, float vwind) {␊
    // Non-BIPV/T walls keep using the legacy wall balance; no special heating term.␊
    if (!wall.hasBipvModel()) {␊
        wall.clearBipvExteriorTemperatures();␊
        wall.recordBipvElectricProduction(0.f);␊
        return 0.f;␊
    }␊
␊
    const WallPVDefinition& pv = *wall.getWallPVDefinition();␊
    const double sigma = 5.670374419e-8;␊
    const double dt = 3600.0; // fixed timestep for the BIPV/T model␊
␊
    // convert temperatures to Kelvin for nonlinear terms␊
    const double TaK = getTa() + 273.15;␊
    const double TwK = Tw + 273.15;␊
    const double ToutK = Tout + 273.15;␊
␊
    const double Awall = max(1e-6, static_cast<double>(wall.getWallArea()));␊
    const double pvratio = static_cast<double>(wall.getPVRatio());␊
    const double L = max(0.1f, wall.getBipvHeight());␊
    const double w = max(0.01f, wall.getBipvWidth());␊
    const double e = max(1e-4, static_cast<double>(pv.eair));␊
␊
    // solar and infrared terms (per wall)␊
    const double Qsun1 = Awall * (wall.getShortWaveIrradiance() * (1.0 - wall.getShortWaveReflectance_opaque())␊
                                  + wall.get_hr() * wall.getEnvironmentalTemperature());␊
    const double envTempK = wall.getEnvironmentalTemperature() + 273.15;␊
    const double wallLongWaveArea = (1.0 - pvratio) * Awall;␊
    const double pvLongWaveArea = pvratio * Awall;␊
    const double hr_env = wall.get_hr();␊
    const double Qsun2 = wall.getGlazingRatio() * wall.getArea()␊
                         * (wall.getGlazingGvalue(wall.getBeamAngle()) * wall.getBeamIrradiance()␊
                            + wall.getGlazingGvalueHemispherical()␊
                              * (wall.getShortWaveIrradiance() - wall.getBeamIrradiance()))␊
                         * (0.15f + 0.85f * wall.getLowerShadingState());␊
␊
    const double rho_air_base = 101325.0 / (287.05 * TaK);␊
␊
    // Newton–Raphson state initialisation (Kelvin for temperatures).␊
    // The convective coefficients start from the façade hc while the airflow begins at the measured␊
    // perpendicular wind or a small positive value to keep Jacobian entries finite. The zone air␊
    // temperature estimate uses the latest Ta() sample so the solver works with current-step values.␊
    array<double,13> x = {ToutK,                             // T_os␊
                          TaK,                               // T_bp␊
                          TaK,                               // T_cell␊
                          TaK,                               // T_gi␊
                          ToutK,                             // T_go␊
                          wall.get_hc(),                     // h_os initial guess␊
                          wall.get_hc(),                     // h_bp initial guess␊
                          TaK,                               // T_air,mean␊
                          TaK,                               // T_air,out␊
                          TaK,                               // X␊
                          0.1,                               // Y␊
                          std::max<double>(0.1, vwind),      // V_air␊
                          0.0};                              // H␊
␊
    auto computeF = [&](const array<double,13>& xin, array<double,13>& F) {␊
        const double Tos = xin[0];␊
        const double Tbp = xin[1];␊
        const double Tcell = xin[2];␊
        const double Tgi = xin[3];␊
        const double Tgo = xin[4];␊
        const double hos = max(0.0, xin[5]);␊
        const double hbp = max(0.0, xin[6]);␊
        const double Tair = xin[7];␊
        const double TairOut = xin[8];␊
        const double Xflow = xin[9];␊
        const double Yflow = xin[10];␊
        const double Vair = max(1e-6, xin[11]);␊
        const double H = xin[12];␊
␊
        // air properties at mean temperature␊
        const double rho_air = 101325.0 / (287.05 * max(1.0, Tair));␊
        const double cp_air = 1005.0;␊
        const double k_air = 2.873e-3 + 7.76e-8 * Tair;␊
        const double nu_air = 3.723e-6 + 4.94e-8 * Tair;␊
        const double beta = 1.0 / max(1.0, Tair);␊
        const double alpha_air = k_air / (rho_air * cp_air);␊
␊
        const double Ra_os = 9.81 * beta * fabs(Tos - Tair) * pow(e, 3.0) / (alpha_air * nu_air);␊
        const double Ra_bp = 9.81 * beta * fabs(Tbp - Tair) * pow(e, 3.0) / (alpha_air * nu_air);␊
␊
        auto natural_hos = [&](double Ra_os) {␊
            if (Ra_os <= 0.) return 0.0;␊
            double ratio_os = Ra_os * e / L;␊
            double invNu = (144.0 / (ratio_os * ratio_os)) + (2.87 / sqrt(max(1e-12, ratio_os)));␊
            return (invNu <= 0.) ? 0.0 : (1.0 / sqrt(invNu)) * k_air / e;␊
        };␊
␊
        auto natural_hbp = [&](double Ra_bp) {␊
            if (Ra_bp <= 0.) return 0.0;␊
            double ratio_bp = Ra_bp * e / L;␊
            double invNu = (144.0 / (ratio_bp * ratio_bp)) + (2.87 / sqrt(max(1e-12, ratio_bp)));␊
            return (invNu <= 0.) ? 0.0 : (1.0 / sqrt(invNu)) * k_air / e;␊
        };␊
␊
␊
        const double hos_nat = natural_hos(Ra_os);␊
        const double hbp_nat = natural_hbp(Ra_bp);␊
␊
        const double Ke = std::max<double>(0.0, getKe());␊
        const double k2 = std::max<double>(1e-6, Ki);␊
        const double Cv = 0.27;␊
        const double Cd = 0.65;␊
␊
        // Equation 1␊
        double QirWall = hr_env * wallLongWaveArea * (envTempK - Tos);␊
        double QirPv = hr_env * pvLongWaveArea * (envTempK - Tgo);␊
␊
        double A11 = pow(Kw1, 2.0) * dt / (max(1e-6f, Cw) + dt * (k2 + Kw1)) - Kw1 - Ke * (1.0 - pvratio);␊
        double B1 = Ke * (1.0 - pvratio) * ToutK + Qsun1 * (1.0 - pvratio) + QirWall
                    + Kw1 / (max(1e-6f, Cw)+dt * (k2 + Kw1)) * (TwK + dt * k2 * TaK + dt * (k2 / max(1e-6, Ki)) * (Qsun2 * Ww + Lr));␊
        double C1 = hos * Awall * pvratio * (Tair - Tos) + pvratio * Awall * sigma * (pv.epsilonbp * pow(Tbp, 4.0) - pv.epsilonos * pow(Tos, 4.0));␊
        F[0] = A11 * Tos + B1 + C1;␊
␊
        // Equation 2␊
        double A22 = -pv.kbp / max(1e-6f, pv.ebp);␊
        double A23 = pv.kbp / max(1e-6f, pv.ebp);␊
        double C2 = hbp * (Tair - Tbp) + sigma * (pv.epsilonos * pow(Tos, 4.0) - pv.epsilonbp * pow(Tbp, 4.0));␊
        F[1] = A22 * Tbp + A23 * Tcell + C2;␊
␊
        // Equation 3␊
        double A32 = (pv.kbp / max(1e-6f, pv.ebp)) * Awall * pvratio;␊
        double A33 = -((pv.alfacsw * pv.taugsw * pv.taugsw * Qsun1 * pvratio + pv.alfacir * pv.taugir * pv.taugir * QirPv) * pv.nuetamp)␊
                     - (pv.kbp / max(1e-6f, pv.ebp)) * Awall * pvratio - (pv.keva / max(1e-6f, pv.eeva)) * Awall * pvratio;␊
        double A34 = (pv.keva / max(1e-6f, pv.eeva)) * Awall * pvratio;␊
        double B3 = (pv.alfacsw * pv.taugsw * pv.taugsw * Qsun1 * pvratio + pv.alfacir * pv.taugir * pv.taugir * QirPv)␊
                    - (pv.alfacsw * pv.taugsw * pv.taugsw * Qsun1 * pvratio + pv.alfacir * pv.taugir * pv.taugir * QirPv)␊
                      * (pv.etaeleref - pv.nuetamp * pv.Tcellref);␊
        F[2] = A32 * Tbp + A33 * Tcell + A34 * Tgi + B3;␊
␊
        // Equation 4␊
        double A43 = pv.keva / max(1e-6f, pv.eeva);␊
        double A44 = -pv.kairg / max(1e-6f, pv.eairg) - pv.keva / max(1e-6f, pv.eeva);␊
        double A45 = pv.kairg / max(1e-6f, pv.eairg);␊
        double C4 = sigma * (pv.epsilongo * pow(Tgo, 4.0) - pv.epsilongi * pow(Tgi, 4.0));␊
        F[3] = A43 * Tcell + A44 * Tgi + A45 * Tgo + C4;␊
␊
        // Equation 5␊
        double A55 = -Ke * pvratio - pv.kairg / max(1e-6f, pv.eairg) * Awall * pvratio;␊
        double A54 = pv.kairg / max(1e-6f, pv.eairg) * Awall * pvratio;␊
        double B5 = Ke * ToutK * pvratio + (pv.algagsw * Qsun1 * pvratio + pv.algagir * QirPv);␊
        double C5 = sigma * Awall * pvratio * (pv.epsilongi * pow(Tgi, 4.0) - pv.epsilongo * pow(Tgo, 4.0));␊
        F[4] = A55 * Tgo + A54 * Tgi + B5 + C5;␊
␊
        // Equations 6 and 7␊
        double A612 = (getTa() <= Tout) ? 0.0 : -0.85 * 1.33;␊
        double A712 = (getTa() <= Tout) ? 0.0 : -0.85 * 1.33;␊
        double B6 = (getTa() <= Tout) ? 0.0 : -0.85 * 1.959;␊
        double B7 = B6;␊
        double C6 = (getTa() <= Tout) ? - hos_nat : -0.85 * 1.517 * pow(fabs(Tos-Tair), 1.0 / 3.0);␊
        double C7 = (getTa() <= Tout) ? - hbp_nat : -0.85 * 1.517 * pow(fabs(Tbp-Tair), 1.0 / 3.0);␊
        F[5] = hos + A612 * Vair + B6 + C6;␊
        F[6] = hbp + A712 * Vair + B7 + C7;␊
␊
        // Equation 8␊
        double C8 = (getTa() <= Tout ? -(TaK - Xflow) : -(ToutK - Xflow)) * ((1.0 - exp(-Yflow * L)) / max(1e-6, Yflow * L));␊
        F[7] = Tair + C8;␊
␊
        // Equation 9␊
        double C9 = (getTa() <= Tout ? -(TaK - Xflow) : -(ToutK - Xflow)) * exp(-Yflow * L);␊
        F[8] = TairOut + C9;␊
␊
        // Equation 10␊
        F[9] = Xflow - (hos * Tos + hbp * Tbp) / max(1e-6, hos + hbp);␊
␊
        // Equation 11␊
        F[10] = Yflow - (hos + hbp) / (e * Vair * rho_air * cp_air);␊
␊
        // Equation 12␊
        double buoyancy = (getTa() <= Tout) ? fabs(Tair - TaK) / max(1e-6, TaK) : fabs(Tair - ToutK) / max(1e-6, ToutK);␊
        double C12 = -Cd * e * w * sqrt(max(0.0, 2.0 * 9.81 * L * buoyancy));␊
        double B12 = (getTa() <= Tout) ? 0.0 : -Cv * e * w * vwind;␊
        F[11] = C12 + B12 + Vair;␊
␊
        // Equation 13 – heating balance for the zone␊
        const double wallDenom = fabs(dt * (-k2 - Kw1) - Cw) < 1e-6 ? (dt * (-k2 - Kw1) - Cw >= 0 ? 1e-6 : -1e-6)␊
                                                                     : (dt * (-k2 - Kw1) - Cw);␊
        double A131 = -(k2 * dt * Kw1) / wallDenom;␊
        const bool heatingMode = (TairOut > TaK) && (TaK < 20.0);␊
        double A1312 = heatingMode ? -e * w * rho_air_base * cp_air * TaK : 0.0;␊
        double B13 = ((-Ci / dt) - getUA() - k2 - (k2 * k2 * dt) / wallDenom)␊
                     * (TaK) + (Ci / dt) * (Ta.back() + 273.15)␊
                     + (getUA() * ToutK + (k2 / std::max<double>(1e-6, Kw1)) * (Qsun2 * Ww + Lr) + Qsun2 * Wa + Lc)
                     - k2 * (Cw * TwK - dt * (k2 / max(1e-6, Ki)) * (Qsun2 * Ww + Lc)) / wallDenom;␊
        double C13 = heatingMode ? e * w * rho_air_base * cp_air * TairOut * Vair : 0.0;␊
        F[12] = A131 * Tos + A1312 * Vair + B13 + C13 + H;␊
    };␊
␊
    array<double,13> F{};␊
␊
    const unsigned int maxIter = 20;␊
    // Newton–Raphson loop with finite-difference Jacobian evaluation.␊
    for (unsigned int iter = 0; iter < maxIter; ++iter) {␊
        computeF(x, F);␊
␊
        bool invalidState = false;␊
        for (double v : F) {␊
            if (!std::isfinite(v)) { invalidState = true; break; }␊
        }␊
        for (double v : x) {␊
            if (!std::isfinite(v)) { invalidState = true; break; }␊
        }␊
        if (invalidState) {␊
            const float nanVal = std::numeric_limits<float>::quiet_NaN();␊
            wall.recordBipvState(nanVal, nanVal, nanVal, nanVal, 0.f);␊
            wall.clearBipvExteriorTemperatures();␊
            wall.recordBipvElectricProduction(0.f);␊
            return 0.f;␊
        }␊
␊
        double norm = 0.0;␊
        for (double v : F) norm = max(norm, fabs(v));␊
        if (norm < 1e-3) break;␊
␊
        double J[13][13];␊
        const double eps = 1e-4;␊
        for (unsigned int j = 0; j < 13; ++j) {␊
            array<double,13> xpert = x;␊
            double delta = eps * max(1.0, fabs(x[j]));␊
            xpert[j] += delta;␊
            array<double,13> Fpert{};␊
            computeF(xpert, Fpert);␊
            for (unsigned int i = 0; i < 13; ++i) {␊
                J[i][j] = (Fpert[i] - F[i]) / delta;␊
            }␊
        }␊
␊
        // solve J * dx = F using Gaussian elimination␊
        double dx[13];␊
        double A[13][14];␊
        for (unsigned int i = 0; i < 13; ++i) {␊
            for (unsigned int j = 0; j < 13; ++j) A[i][j] = J[i][j];␊
            A[i][13] = F[i];␊
        }␊
␊
        for (unsigned int k = 0; k < 13; ++k) {␊
            // pivot␊
            unsigned int pivot = k;␊
            for (unsigned int i = k + 1; i < 13; ++i) {␊
                if (fabs(A[i][k]) > fabs(A[pivot][k])) pivot = i;␊
            }␊
            if (fabs(A[pivot][k]) < 1e-12) continue;␊
            if (pivot != k) {␊
                for (unsigned int j = k; j < 14; ++j) swap(A[k][j], A[pivot][j]);␊
            }␊
            double diag = A[k][k];␊
            for (unsigned int j = k; j < 14; ++j) A[k][j] /= diag;␊
            for (unsigned int i = 0; i < 13; ++i) {␊
                if (i == k) continue;␊
                double factor = A[i][k];␊
                for (unsigned int j = k; j < 14; ++j) A[i][j] -= factor * A[k][j];␊
            }␊
        }␊
␊
        for (unsigned int i = 0; i < 13; ++i) dx[i] = A[i][13];␊
␊
        for (unsigned int i = 0; i < 13; ++i) x[i] -= dx[i];␊
    }␊
␊
    wall.setBipvExteriorTemperatures(x[0], x[4]);␊
␊
    const double QirPvTerm = hr_env * pvLongWaveArea * (envTempK - x[4]);␊
␊
    // BIPV/T electrical production (per timestep) follows the explicit balance:␊
    //   E_el = pvratio * (α_c,sw τ_g,sw^2 Q_sun1 + α_c,ir τ_g,ir^2 Q_ir,PV) * η_el,␊
    //   η_el = max(0, η_el,ref − ν_etamp (T_cell − T_cell,ref)).␊
    // The first term represents the irradiance incident on the PV fraction of the façade␊
    // (shortwave + longwave), while η_el applies the temperature-dependent electrical efficiency.␊
    const double incidentSolar = pvratio * (pv.alfacsw * pv.taugsw * pv.taugsw * Qsun1␊
                                            + pv.alfacir * pv.taugir * pv.taugir * QirPvTerm);␊
    double etaElectric = pv.etaeleref - pv.nuetamp * (x[2] - pv.Tcellref);␊
    etaElectric = std::max(0.0, etaElectric);␊
    wall.recordBipvState(static_cast<float>(x[0]), static_cast<float>(x[4]), static_cast<float>(x[2]), static_cast<float>(x[1]), static_cast<float>(etaElectric));␊
    wall.recordBipvElectricProduction(static_cast<float>(incidentSolar * etaElectric));␊
␊
    return static_cast<float>(x[12]);␊
}␊
␊
// Zone3N␊
Zone3N::Zone3N(unsigned int id, Building* pBuilding, bool groundFloor, float Vi, vector<Wall*> walls, vector<Roof*> roofs, vector<Surface*> surfaces, vector<Floor*> floors, Occupants *pOccupants)␊
        :Zone2N(id,pBuilding,groundFloor,Vi,walls,roofs,surfaces,floors,pOccupants) {␊
␊
    // intialisation␊
    nNodes=3;␊
    Tr=15.f;␊
␊
    update(true);␊
}␊
␊
/**␊
 * @brief Zone3N::updateSimulationModelParameters: update all "secondary" parameters, i.e.␊
 *  parameters which are based on other parameters that can be modified in the interface CitySimPro␊
 */␊
void Zone3N::update(bool constructor){␊
    if(!constructor){␊
        Zone2N::update();␊
    }␊
    // If constructor, this was already called when the parent Zone2N was initialised␊
␊
    if (roofs.empty()) { // no roofs in this zone, temperature decoupled from the model␊
        Cr = 1.f;␊
        Kr1 = 1.f;␊
        Kr2 = 1.f;␊
    }␊
    else {␊
        roofs[0]->getComposite()->getSimplifiedNode(Cr, Kr1, Kr2);␊
        // multiply with the total roofs area␊
        Cr *= Sro;␊
    }␊
␊
    // outputs the values␊
    //logStream << "Cr: " << Cr << "\tKr1: " << Kr1*Sro << "\tKr2: " << Kr2*Sro << endl;␊
}␊
␊
double Zone3N::getKappa3() {␊
␊
    // if no roofs, returns a dummy value to couple with the outdoor environment (avoid singular matrix)␊
    if (roofs.empty()) return 1.;␊
␊
    double kappa3 = 0.;␊
    for (size_t i=0; i<roofs.size(); ++i){␊
        kappa3 += Kr1*roofs[i]->getRoofArea()*(roofs[i]->get_hc()+roofs[i]->get_hr()+roofs[i]->get_X())/(Kr1+roofs[i]->get_hc()+roofs[i]->get_hr()+roofs[i]->get_X());␊
        //cout << "Zone3N::getKappa3: hc="<< roofs[i]->get_hc()<<" hr="<<roofs[i]->get_hr() << " area="<<roofs[i]->getRoofArea()<<endl<<flush;␊
    }␊
    return kappa3;␊
␊
}␊
␊
void Zone3N::setTos(float Tout) {␊
    // define the surface temperature per surface␊
    for (size_t i=0; i<walls.size(); ++i) {␊
        float wallTemperature = ( Kw1*Tw + walls[i]->get_hc()*Tout␊
                                  +walls[i]->getShortWaveIrradiance()*(1.-walls[i]->getShortWaveReflectance_opaque())␊
                                  +walls[i]->get_hr()*walls[i]->getEnvironmentalTemperature() )␊
                                / (Kw1 + walls[i]->get_hc() + walls[i]->get_hr());␊
        if (wallTemperature>100) {␊
            logStream << "Wall: " + toString(walls[i]->getId()) + "(" + toString(walls[i]->getKey()) + "), Temperature: " + toString(wallTemperature)␊
                      << "\tTw: " << Tw << "\tTout: " << Tout << "\tIrradiance: " << walls[i]->getShortWaveIrradiance() << "\tIrradiance absorbed: " << walls[i]->getShortWaveIrradiance()*(1.-walls[i]->getShortWaveReflectance_opaque())␊
                      << "\tEnvironmental Temperature: " << walls[i]->get_hr()*walls[i]->getEnvironmentalTemperature() << endl;␊
        }␊
        walls[i]->setTemperature(wallTemperature);␊
    }␊
    for (size_t i=0; i<roofs.size(); ++i) {␊
        float roofTemperature = ( Kr1*Tr + roofs[i]->get_hc()*Tout␊
                                  +roofs[i]->getShortWaveIrradiance()*(1.-roofs[i]->getShortWaveReflectance_opaque())␊
                                  +roofs[i]->get_hr()*roofs[i]->getEnvironmentalTemperature()␊
                                  -roofs[i]->get_Y() )␊
                                / (Kr1 + roofs[i]->get_hc() + roofs[i]->get_hr() + roofs[i]->get_X());␊
        if (roofTemperature>100) {␊
            logStream << "Roof: " + toString(roofs[i]->getId()) + "(" + toString(roofs[i]->getKey()) + "), Temperature: " + toString(roofTemperature)␊
                      << "Tr: " << Tr << "\tTout: " << Tout << "\tIrradiance: " << roofs[i]->getShortWaveIrradiance() << "\tIrradiance absorbed: " << roofs[i]->getShortWaveIrradiance()*(1.-roofs[i]->getShortWaveReflectance_opaque())␊
                      << "\tEnvironmental Temperature: " << roofs[i]->get_hr()*roofs[i]->getEnvironmentalTemperature() << endl;␊
        }␊
        roofs[i]->setTemperature(roofTemperature);␊
    }␊
}␊
␊
// Zone3N␊
Zone3N_floor::Zone3N_floor(unsigned int id, Building* pBuilding, bool groundFloor, float Vi, vector<Wall*> walls, vector<Roof*> roofs, vector<Surface*> surfaces, vector<Floor*> floors, Occupants *pOccupants)␊
        :Zone2N(id,pBuilding,groundFloor,Vi,walls,roofs,surfaces,floors,pOccupants) {␊
␊
    // intialisation␊
    nNodes=3;␊
    Tf=15.f;␊
␊
    update(true);␊
␊
}␊
␊
/**␊
 * @brief Zone3N_floor::updateSimulationModelParameters: update all "secondary" parameters, i.e.␊
 *  parameters which are based on other parameters that can be modified in the interface CitySimPro␊
 */␊
void Zone3N_floor::update(bool constructor){␊
    if(!constructor){␊
        Zone2N::update();␊
    }␊
    // If constructor, this was already called when the parent Zone2N was initialised␊
␊
    if (floors.empty()) { // no floor in this zone␊
        throw(string("Creation of a thermal zone without floor using the four nodes model."));␊
        //Cw = 1.;␊
        //Kw1 = 1.;␊
        //Kw2 = 1.;␊
        //Ki = 0.; // the computation of the wall temperature is completely disconnected from the air node temperature␊
    }␊
    else floors[0]->getComposite()->getSimplifiedNode(Cf, Kf1, Kf2);␊
␊
    // multiply with the total floor area␊
    Sf = 0.f;␊
    for (size_t i=0; i<floors.size(); ++i) Sf += floors[i]->getArea();␊
    Cf *= Sf;␊
␊
    // outputs the values␊
    //logStream << "Cf: " << Cf << "\tKf1: " << Kf1*Sf << "\tKf2: " << Kf2*Sf << endl;␊
}␊
␊
// Zone4N␊
Zone4N::Zone4N(unsigned int id, Building* pBuilding, bool groundFloor, float Vi, vector<Wall*> walls, vector<Roof*> roofs, vector<Surface*> surfaces, vector<Floor*> floors, Occupants *pOccupants)␊
        :Zone3N(id,pBuilding,groundFloor,Vi,walls,roofs,surfaces,floors,pOccupants) {␊
␊
    logStream << "Zone4N with 4 thermal nodes." << endl;␊
␊
    // intialisation␊
    nNodes=4;␊
    Tf=15.f;␊
␊
    update(true);␊
}␊
␊
Zone4N::Zone4N(unsigned int id, Building* pBuilding, bool groundFloor, float Vi, Surface* floor, float elevation)␊
    :Zone3N(id,pBuilding,groundFloor,Vi,vector<Wall*>(),vector<Roof*>(),vector<Surface*>(),vector<Floor*>(),NULL){␊
    // intialisation␊
    nNodes=4;␊
    Tf=15.f;␊
␊
    floors.push_back(new Floor(*floor));␊
    roofs.push_back(new Roof(*floor, elevation)); // will get floor id + 1␊
    vector<GENPoint>* vertices = floor->getVertices();␊
    for(unsigned int i=0; i<vertices->size()-1; ++i){␊
        walls.push_back(new Wall(floor->getId()+1+i,(*vertices)[i],(*vertices)[i+1],elevation));␊
    }␊
    walls.push_back(new Wall(floor->getId()+vertices->size(),(*vertices)[vertices->size()-1],(*vertices)[0],elevation));␊
    delete floor;␊
}␊
␊
Zone4N::Zone4N(Building* pBuilding, bool groundFloor):Zone3N(pBuilding->getId(),pBuilding,groundFloor){␊
    // Incomplete constructor for DXF reading, do not use without completing the geometry...␊
    nNodes=4;␊
    Tf=15.f;␊
}␊
␊
void Zone4N::addSurface(Surface* s){␊
    //cout << "Zone::addSurface, normal altitude: " << s->normal().Altitude().degrees() << endl;␊
    if(s->normal().Altitude().degrees() > 1.){ // roof␊
        //cout << "Zone::addSurface -> roof" << endl;␊
        roofs.push_back(new Roof(*s));␊
    }␊
    else if(s->normal().Altitude().degrees() < -1.){ // floor␊
        //cout << "Zone::addSurface -> floor" << endl;␊
        floors.push_back(new Floor(*s));␊
    }␊
    else{ // all other surfaces added through this function are considered walls␊
        //cout << "Zone::addSurface -> wall" << endl;␊
        walls.push_back(new Wall(*s));␊
    }␊
    delete s;␊
}␊
␊
/**␊
 * @brief Zone4N::updateSimulationModelParameters: update all "secondary" parameters, i.e.␊
 *  parameters which are based on other parameters that can be modified in the interface CitySimPro␊
 */␊
void Zone4N::update(bool constructor){␊
    if(!constructor){␊
        Zone3N::update();␊
    }␊
    // If constructor, this was already called when the parent Zone3N was initialised␊
␊
    if (floors.empty()) { // no floor in this zone, temperature decoupled from the model␊
        Cf = 1.f;␊
        Kf1 = 1.f;␊
        Kf2 = 1.f;␊
    }␊
    else  {␊
        floors[0]->getComposite()->getSimplifiedNode(Cf, Kf1, Kf2);␊
        // multiply with the total floor area␊
        Sf = 0.f;␊
        for (size_t i=0; i<floors.size(); ++i) Sf += floors[i]->getArea();␊
        Cf *= Sf;␊
    }␊
␊
    // outputs the values␊
    //logStream << "Cf: " << Cf << "\tKf1: " << Kf1*Sf << "\tKf2: " << Kf2*Sf << endl;␊
}␊
␊
double Zone4N::getKappa5() {␊
␊
    // if no floors, returns a dummy value to couple with the outdoor environment (avoid singular matrix)␊
    if (floors.empty()) return 1.;␊
    else return Kf1*Sf;␊
␊
}␊
␊
// ZoneN␊
ZoneN::ZoneN(unsigned int id, Building* pBuilding, bool groundFloor, float Vi, vector<Wall*> walls, vector<Roof*> roofs, vector<Surface*> surfaces, vector<Floor*> floors, Occupants *pOccupants)␊
        :Zone(id,pBuilding,groundFloor,Vi,walls,roofs,surfaces,floors,pOccupants) {␊
␊
    update(true);␊
}␊
␊
/**␊
 * @brief ZoneN::updateSimulationModelParameters: update all "secondary" parameters, i.e.␊
 *  parameters which are based on other parameters that can be modified in the interface CitySimPro␊
 */␊
void ZoneN::update(bool constructor){␊
    if(!constructor){␊
        Zone::update();␊
    }␊
    // If constructor, this was already called when the parent ZoneN was initialised␊
␊
    // check the potential errors␊
    if (walls.empty()) { // no external wall in this zone␊
        throw(string("Creation of a thermal zone without external walls using the n-nodes model."));␊
        //Cw = 1.;␊
        //Kw1 = 1.;␊
        //Kw2 = 1.;␊
        //Ki = 0.; // the computation of the wall temperature is completely disconnected from the air node temperature␊
    }␊
    else {␊
        // intialisation␊
        nNodes=walls[0]->getComposite()->getnLayers()+1; // plus one node for the internal air␊
    }␊
␊
    Tw.assign(nNodes-1, 15.f);␊
    if (Model::thermalExplicit) TwExpl.assign(nNodes-1, 15.f);␊
␊
    // initialisation of the Ki (using a fixed value of hc=3.0)␊
    Ki=3.0*Swa;␊
}␊
␊
void ZoneN::setTos(float Tout) {␊
    // clamp the wall temperature, at minimum the outside air temperature, at maximum 80∞C (upper limit value taken from CIBSE guide)␊
    //float wallTemperature = min(max((Kw1*(Tw.back()) + getKe()*(Tout) + getQsun1())/(Kw1 + getKe()), Tout), 80.);␊
    float wallTemperature = (getG0()*Tw[0] + getKe()*Tout + getQsun1())/(getG0() + getKe() + getHr());␊
    for (size_t i=0; i<walls.size(); ++i) {␊
        // sets the same wall temperature for all walls as only one temperature for the walls␊
        walls[i]->setTemperature(wallTemperature);␊
    }␊
    for (size_t i=0; i<roofs.size(); ++i) {␊
        // clamp the roof temperature, at minimum the outside air temperature, at maximum 80∞C (upper limit value taken from CIBSE guide)␊
        //float roofTemperature = min(max(( Kr*Ta.back() + (getKe()/getSwa())*Tout + getQsun3(i)/roofs[i]->getRoofArea())/(Kr + 23.f), Tout), 80.);␊
        float roofTemperature = ( roofs[i]->getKr()*Ta.back() + roofs[i]->get_hc()*Tout␊
                                  +roofs[i]->getShortWaveIrradiance()*(1.-roofs[i]->getShortWaveReflectance_opaque())␊
                                  +roofs[i]->get_hr()*roofs[i]->getEnvironmentalTemperature() )␊
                                / (roofs[i]->getKr() + roofs[i]->get_hc() + roofs[i]->get_hr());␊
        roofs[i]->setTemperature(roofTemperature);␊
    }␊
}␊
␊
␊
␊

