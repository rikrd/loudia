/*                                                         
** Copyright (C) 2008, 2009 Ricard Marxer <email@ricardmarxer.com>
**                                                                  
** This program is free software; you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation; either version 2 of the License, or   
** (at your option) any later version.                                 
**                                                                     
** This program is distributed in the hope that it will be useful,     
** but WITHOUT ANY WARRANTY; without even the implied warranty of      
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       
** GNU General Public License for more details.                        
**                                                                     
** You should have received a copy of the GNU General Public License   
** along with this program; if not, write to the Free Software         
** Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307, USA.
*/

#include "FilterUtils.h"

using namespace std;

void chebyshev1(int order, Real rippleDB, int channels, MatrixXC* zeros, MatrixXC* poles, Real* gain) {
  (*zeros) = MatrixXC::Zero(channels, 1);
 
  Real eps = sqrt(pow(10, (0.1 * rippleDB)) - 1.0);

  MatrixXC n;
  range(1, order + 1, order, channels, &n);

  Real mu = 1.0 / order * log((1.0 + sqrt( 1 + eps * eps)) / eps);

  MatrixXC theta = ((n * 2).cwise() - 1.0) / order * M_PI / 2.0;

  (*poles) = -sinh(mu) * theta.cwise().sin() + Complex(0, 1) * cosh(mu) * theta.cwise().cos();

  Complex gainComplex = 1.0;

  for ( int i = 0; i < (*poles).cols(); i++ ) {
    gainComplex *= -(*poles)(0, i);
  }

  (*gain) = gainComplex.real();
  
  if ( order % 2 == 0 ) {
    (*gain) /= sqrt((1 + eps * eps));
  }
}

void chebyshev2(int order, Real rippleDB, int channels, MatrixXC* zeros, MatrixXC* poles, Real* gain) {
  Real de = 1.0 / sqrt(pow(10, (0.1 * rippleDB)) - 1.0);
  Real mu = asinh(1.0 / de) / order;

  MatrixXC n;
  if(order % 2) {
    n.resize(channels, order - 1);
    MatrixXC nFirst;
    range(1, order , order/2, channels, &nFirst);

    MatrixXC nSecond;
    range(order + 2, (2 * order) + 1, order/2, channels, &nSecond);
    
    n << nFirst, nSecond;
  } else{
    n.resize(channels, order);
    range(1, (2 * order) + 1, order, channels, &n);
  }
  
  (*zeros) = (Complex(0,1) * ((n * M_PI) / (2.0 * order)).cwise().cos().cwise().inverse()).conjugate();

  MatrixXC rng;
  range(1, (2 * order) + 1, order, channels, &rng);
  
  (*poles) = (Complex(0,1) * (((M_PI * rng) / (2.0*order)).cwise() + M_PI / 2.0)).cwise().exp();

  (*poles) = (((*poles).real().cast<Complex>() * sinh( mu )) + (Complex(0, 1) * cosh( mu ) * (*poles).imag().cast<Complex>())).cwise().inverse();

  // TODO: gain should be a vector (one gain per channel)
  (*gain) = ((-(*poles)).rowwise().prod().cwise() / (-(*zeros)).rowwise().prod()).real().sum();
}

void butterworth(int order, int channels, MatrixXC* zeros, MatrixXC* poles, Real* gain) {
  MatrixXC n;
  range(1, order + 1, order + 1, &n);

  (*zeros) = MatrixXC::Zero(channels, 1);
  (*poles) = (((2*n).cwise() - 1) * Complex(0, 1) / (2.0 * order) * M_PI).cwise().exp() * Complex(0, 1);
  (*gain) = 1.0;
}

void bessel(int order, int channels, MatrixXC* zeros, MatrixXC* poles, Real* gain) {
  (*gain) = 1;
  (*zeros) = MatrixXC::Zero(channels, 1);
  MatrixXC tempPoles(1, order);
  
  switch( order ) {
  case 1:
    tempPoles << -1;
    break;
    
  case 2:
    tempPoles <<  Complex(-.8660254037844386467637229, .4999999999999999999999996),
      Complex(-.8660254037844386467637229, -.4999999999999999999999996);
    break;

  case 3:
    tempPoles << Complex(-.9416000265332067855971980, 0),
      Complex(-.7456403858480766441810907, -.7113666249728352680992154),
      Complex(-.7456403858480766441810907, +.7113666249728352680992154);
    break;

  case 4:
    tempPoles << Complex(-.6572111716718829545787781, -.8301614350048733772399715),
      Complex(-.6572111716718829545787788, +.8301614350048733772399715),
      Complex(-.9047587967882449459642637, -.2709187330038746636700923),
      Complex(-.9047587967882449459642624, +.2709187330038746636700926);
    break;

  case 5:
    tempPoles << Complex(-.9264420773877602247196260, 0),
      Complex(-.8515536193688395541722677, -.4427174639443327209850002),
      Complex(-.8515536193688395541722677, +.4427174639443327209850002),
      Complex(-.5905759446119191779319432, -.9072067564574549539291747),
      Complex(-.5905759446119191779319432, +.9072067564574549539291747);
    break;
    
  case 6:
    tempPoles << Complex(-.9093906830472271808050953, -.1856964396793046769246397),
      Complex(-.9093906830472271808050953, +.1856964396793046769246397),
      Complex(-.7996541858328288520243325, -.5621717346937317988594118),
      Complex(-.7996541858328288520243325, +.5621717346937317988594118),
      Complex(-.5385526816693109683073792, -.9616876881954277199245657),
      Complex(-.5385526816693109683073792, +.9616876881954277199245657);
    break;

  case 7:
    tempPoles << Complex(-.9194871556490290014311619, 0),
      Complex(-.8800029341523374639772340, -.3216652762307739398381830),
      Complex(-.8800029341523374639772340, +.3216652762307739398381830),
      Complex(-.7527355434093214462291616, -.6504696305522550699212995),
      Complex(-.7527355434093214462291616, +.6504696305522550699212995),
      Complex(-.4966917256672316755024763, -1.002508508454420401230220),
      Complex(-.4966917256672316755024763, +1.002508508454420401230220);
    break;

  case 8:
    tempPoles << Complex(-.9096831546652910216327629-.1412437976671422927888150),
      Complex(-.9096831546652910216327629, +.1412437976671422927888150),
      Complex(-.8473250802359334320103023, -.4259017538272934994996429),
      Complex(-.8473250802359334320103023, +.4259017538272934994996429),
      Complex(-.7111381808485399250796172, -.7186517314108401705762571),
      Complex(-.7111381808485399250796172, +.7186517314108401705762571),
      Complex(-.4621740412532122027072175, -1.034388681126901058116589),
      Complex(-.4621740412532122027072175, +1.034388681126901058116589);
    break;
    
  case 9:
    tempPoles << Complex(-.9154957797499037686769223, 0),
      Complex(-.8911217017079759323183848, -.2526580934582164192308115),
      Complex(-.8911217017079759323183848, +.2526580934582164192308115),
      Complex(-.8148021112269012975514135, -.5085815689631499483745341),
      Complex(-.8148021112269012975514135, +.5085815689631499483745341),
      Complex(-.6743622686854761980403401, -.7730546212691183706919682),
      Complex(-.6743622686854761980403401, +.7730546212691183706919682),
      Complex(-.4331415561553618854685942, -1.060073670135929666774323),
      Complex(-.4331415561553618854685942, +1.060073670135929666774323);
    break;
    
  case 10:
    tempPoles << Complex(-.9091347320900502436826431, -.1139583137335511169927714),
      Complex(-.9091347320900502436826431, +.1139583137335511169927714),
      Complex(-.8688459641284764527921864, -.3430008233766309973110589),
      Complex(-.8688459641284764527921864, +.3430008233766309973110589),
      Complex(-.7837694413101441082655890, -.5759147538499947070009852),
      Complex(-.7837694413101441082655890, +.5759147538499947070009852),
      Complex(-.6417513866988316136190854, -.8175836167191017226233947),
      Complex(-.6417513866988316136190854, +.8175836167191017226233947),
      Complex(-.4083220732868861566219785, -1.081274842819124562037210),
      Complex(-.4083220732868861566219785, +1.081274842819124562037210);
    break;
    
  case 11:
    tempPoles << Complex(-.9129067244518981934637318, 0),
      Complex(-.8963656705721166099815744, -.2080480375071031919692341),
      Complex(-.8963656705721166099815744, +.2080480375071031919692341),
      Complex(-.8453044014712962954184557, -.4178696917801248292797448),
      Complex(-.8453044014712962954184557, +.4178696917801248292797448),
      Complex(-.7546938934722303128102142, -.6319150050721846494520941),
      Complex(-.7546938934722303128102142, +.6319150050721846494520941),
      Complex(-.6126871554915194054182909, -.8547813893314764631518509),
      Complex(-.6126871554915194054182909, +.8547813893314764631518509),
      Complex(-.3868149510055090879155425, -1.099117466763120928733632),
      Complex(-.3868149510055090879155425, +1.099117466763120928733632);
    break;
    
  case 12:
    tempPoles << Complex(-.9084478234140682638817772, -95506365213450398415258360.0e-27),
      Complex(-.9084478234140682638817772, +95506365213450398415258360.0e-27),
      Complex(-.8802534342016826507901575, -.2871779503524226723615457),
      Complex(-.8802534342016826507901575, +.2871779503524226723615457),
      Complex(-.8217296939939077285792834, -.4810212115100676440620548),
      Complex(-.8217296939939077285792834, +.4810212115100676440620548),
      Complex(-.7276681615395159454547013, -.6792961178764694160048987),
      Complex(-.7276681615395159454547013, +.6792961178764694160048987),
      Complex(-.5866369321861477207528215, -.8863772751320727026622149),
      Complex(-.5866369321861477207528215, +.8863772751320727026622149),
      Complex(-.3679640085526312839425808, -1.114373575641546257595657),
      Complex(-.3679640085526312839425808, +1.114373575641546257595657);
    break;
    
  case 13:
    tempPoles << Complex(-.9110914665984182781070663, 0),
      Complex(-.8991314665475196220910718, -.1768342956161043620980863),
      Complex(-.8991314665475196220910718, +.1768342956161043620980863),
      Complex(-.8625094198260548711573628, -.3547413731172988997754038),
      Complex(-.8625094198260548711573628, +.3547413731172988997754038),
      Complex(-.7987460692470972510394686, -.5350752120696801938272504),
      Complex(-.7987460692470972510394686, +.5350752120696801938272504),
      Complex(-.7026234675721275653944062, -.7199611890171304131266374),
      Complex(-.7026234675721275653944062, +.7199611890171304131266374),
      Complex(-.5631559842430199266325818, -.9135900338325109684927731),
      Complex(-.5631559842430199266325818, +.9135900338325109684927731),
      Complex(-.3512792323389821669401925, -1.127591548317705678613239),
      Complex(-.3512792323389821669401925, +1.127591548317705678613239);
    break;
    
  case 14:
    tempPoles << Complex(-.9077932138396487614720659, -82196399419401501888968130.0e-27),
      Complex(-.9077932138396487614720659, +82196399419401501888968130.0e-27),
      Complex(-.8869506674916445312089167, -.2470079178765333183201435),
      Complex(-.8869506674916445312089167, +.2470079178765333183201435),
      Complex(-.8441199160909851197897667, -.4131653825102692595237260),
      Complex(-.8441199160909851197897667, +.4131653825102692595237260),
      Complex(-.7766591387063623897344648, -.5819170677377608590492434),
      Complex(-.7766591387063623897344648, +.5819170677377608590492434),
      Complex(-.6794256425119233117869491, -.7552857305042033418417492),
      Complex(-.6794256425119233117869491, +.7552857305042033418417492),
      Complex(-.5418766775112297376541293, -.9373043683516919569183099),
      Complex(-.5418766775112297376541293, +.9373043683516919569183099),
      Complex(-.3363868224902037330610040, -1.139172297839859991370924),
      Complex(-.3363868224902037330610040, +1.139172297839859991370924);
    break;
    
  case 15:
    tempPoles << Complex(-.9097482363849064167228581, 0),
      Complex(-.9006981694176978324932918, -.1537681197278439351298882),
      Complex(-.9006981694176978324932918, +.1537681197278439351298882),
      Complex(-.8731264620834984978337843, -.3082352470564267657715883),
      Complex(-.8731264620834984978337843, +.3082352470564267657715883),
      Complex(-.8256631452587146506294553, -.4642348752734325631275134),
      Complex(-.8256631452587146506294553, +.4642348752734325631275134),
      Complex(-.7556027168970728127850416, -.6229396358758267198938604),
      Complex(-.7556027168970728127850416, +.6229396358758267198938604),
      Complex(-.6579196593110998676999362, -.7862895503722515897065645),
      Complex(-.6579196593110998676999362, +.7862895503722515897065645),
      Complex(-.5224954069658330616875186, -.9581787261092526478889345),
      Complex(-.5224954069658330616875186, +.9581787261092526478889345),
      Complex(-.3229963059766444287113517, -1.149416154583629539665297),
      Complex(-.3229963059766444287113517, +1.149416154583629539665297);
    break;
    
  case 16:
    tempPoles << Complex(-.9072099595087001356491337, -72142113041117326028823950.0e-27),
      Complex(-.9072099595087001356491337, +72142113041117326028823950.0e-27),
      Complex(-.8911723070323647674780132, -.2167089659900576449410059),
      Complex(-.8911723070323647674780132, +.2167089659900576449410059),
      Complex(-.8584264231521330481755780, -.3621697271802065647661080),
      Complex(-.8584264231521330481755780, +.3621697271802065647661080),
      Complex(-.8074790293236003885306146, -.5092933751171800179676218),
      Complex(-.8074790293236003885306146, +.5092933751171800179676218),
      Complex(-.7356166304713115980927279, -.6591950877860393745845254),
      Complex(-.7356166304713115980927279, +.6591950877860393745845254),
      Complex(-.6379502514039066715773828, -.8137453537108761895522580),
      Complex(-.6379502514039066715773828, +.8137453537108761895522580),
      Complex(-.5047606444424766743309967, -.9767137477799090692947061),
      Complex(-.5047606444424766743309967, +.9767137477799090692947061),
      Complex(-.3108782755645387813283867, -1.158552841199330479412225),
      Complex(-.3108782755645387813283867, +1.158552841199330479412225);
    break;
    
  case 17:
    tempPoles << Complex(-.9087141161336397432860029, 0),
      Complex(-.9016273850787285964692844, -.1360267995173024591237303),
      Complex(-.9016273850787285964692844, +.1360267995173024591237303),
      Complex(-.8801100704438627158492165, -.2725347156478803885651973),
      Complex(-.8801100704438627158492165, +.2725347156478803885651973),
      Complex(-.8433414495836129204455491, -.4100759282910021624185986),
      Complex(-.8433414495836129204455491, +.4100759282910021624185986),
      Complex(-.7897644147799708220288138, -.5493724405281088674296232),
      Complex(-.7897644147799708220288138, +.5493724405281088674296232),
      Complex(-.7166893842372349049842743, -.6914936286393609433305754),
      Complex(-.7166893842372349049842743, +.6914936286393609433305754),
      Complex(-.6193710717342144521602448, -.8382497252826992979368621),
      Complex(-.6193710717342144521602448, +.8382497252826992979368621),
      Complex(-.4884629337672704194973683, -.9932971956316781632345466),
      Complex(-.4884629337672704194973683, +.9932971956316781632345466),
      Complex(-.2998489459990082015466971, -1.166761272925668786676672),
      Complex(-.2998489459990082015466971, +1.166761272925668786676672);
    break;
    
  case 18:
    tempPoles << Complex(-.9067004324162775554189031, -64279241063930693839360680.0e-27),
      Complex(-.9067004324162775554189031, +64279241063930693839360680.0e-27),
      Complex(-.8939764278132455733032155, -.1930374640894758606940586),
      Complex(-.8939764278132455733032155, +.1930374640894758606940586),
      Complex(-.8681095503628830078317207, -.3224204925163257604931634),
      Complex(-.8681095503628830078317207, +.3224204925163257604931634),
      Complex(-.8281885016242836608829018, -.4529385697815916950149364),
      Complex(-.8281885016242836608829018, +.4529385697815916950149364),
      Complex(-.7726285030739558780127746, -.5852778162086640620016316),
      Complex(-.7726285030739558780127746, +.5852778162086640620016316),
      Complex(-.6987821445005273020051878, -.7204696509726630531663123),
      Complex(-.6987821445005273020051878, +.7204696509726630531663123),
      Complex(-.6020482668090644386627299, -.8602708961893664447167418),
      Complex(-.6020482668090644386627299, +.8602708961893664447167418),
      Complex(-.4734268069916151511140032, -1.008234300314801077034158),
      Complex(-.4734268069916151511140032, +1.008234300314801077034158),
      Complex(-.2897592029880489845789953, -1.174183010600059128532230),
      Complex(-.2897592029880489845789953, +1.174183010600059128532230);
    break;
    
  case 19:
    tempPoles << Complex(-.9078934217899404528985092, 0),
      Complex(-.9021937639390660668922536, -.1219568381872026517578164),
      Complex(-.9021937639390660668922536, +.1219568381872026517578164),
      Complex(-.8849290585034385274001112, -.2442590757549818229026280),
      Complex(-.8849290585034385274001112, +.2442590757549818229026280),
      Complex(-.8555768765618421591093993, -.3672925896399872304734923),
      Complex(-.8555768765618421591093993, +.3672925896399872304734923),
      Complex(-.8131725551578197705476160, -.4915365035562459055630005),
      Complex(-.8131725551578197705476160, +.4915365035562459055630005),
      Complex(-.7561260971541629355231897, -.6176483917970178919174173),
      Complex(-.7561260971541629355231897, +.6176483917970178919174173),
      Complex(-.6818424412912442033411634, -.7466272357947761283262338),
      Complex(-.6818424412912442033411634, +.7466272357947761283262338),
      Complex(-.5858613321217832644813602, -.8801817131014566284786759),
      Complex(-.5858613321217832644813602, +.8801817131014566284786759),
      Complex(-.4595043449730988600785456, -1.021768776912671221830298),
      Complex(-.4595043449730988600785456, +1.021768776912671221830298),
      Complex(-.2804866851439370027628724, -1.180931628453291873626003),
      Complex(-.2804866851439370027628724, +1.180931628453291873626003);
    break;
    
  case 20:
    tempPoles << Complex(-.9062570115576771146523497, -57961780277849516990208850.0e-27),
      Complex(-.9062570115576771146523497, +57961780277849516990208850.0e-27),
      Complex(-.8959150941925768608568248, -.1740317175918705058595844),
      Complex(-.8959150941925768608568248, +.1740317175918705058595844),
      Complex(-.8749560316673332850673214, -.2905559296567908031706902),
      Complex(-.8749560316673332850673214, +.2905559296567908031706902),
      Complex(-.8427907479956670633544106, -.4078917326291934082132821),
      Complex(-.8427907479956670633544106, +.4078917326291934082132821),
      Complex(-.7984251191290606875799876, -.5264942388817132427317659),
      Complex(-.7984251191290606875799876, +.5264942388817132427317659),
      Complex(-.7402780309646768991232610, -.6469975237605228320268752),
      Complex(-.7402780309646768991232610, +.6469975237605228320268752),
      Complex(-.6658120544829934193890626, -.7703721701100763015154510),
      Complex(-.6658120544829934193890626, +.7703721701100763015154510),
      Complex(-.5707026806915714094398061, -.8982829066468255593407161),
      Complex(-.5707026806915714094398061, +.8982829066468255593407161),
      Complex(-.4465700698205149555701841, -1.034097702560842962315411),
      Complex(-.4465700698205149555701841, +1.034097702560842962315411),
      Complex(-.2719299580251652601727704, -1.187099379810885886139638),
      Complex(-.2719299580251652601727704, +1.187099379810885886139638);
    break;
    
  case 21:
    tempPoles << Complex(-.9072262653142957028884077, 0),
      Complex(-.9025428073192696303995083, -.1105252572789856480992275),
      Complex(-.9025428073192696303995083, +.1105252572789856480992275),
      Complex(-.8883808106664449854431605, -.2213069215084350419975358),
      Complex(-.8883808106664449854431605, +.2213069215084350419975358),
      Complex(-.8643915813643204553970169, -.3326258512522187083009453),
      Complex(-.8643915813643204553970169, +.3326258512522187083009453),
      Complex(-.8299435470674444100273463, -.4448177739407956609694059),
      Complex(-.8299435470674444100273463, +.4448177739407956609694059),
      Complex(-.7840287980408341576100581, -.5583186348022854707564856),
      Complex(-.7840287980408341576100581, +.5583186348022854707564856),
      Complex(-.7250839687106612822281339, -.6737426063024382240549898),
      Complex(-.7250839687106612822281339, +.6737426063024382240549898),
      Complex(-.6506315378609463397807996, -.7920349342629491368548074),
      Complex(-.6506315378609463397807996, +.7920349342629491368548074),
      Complex(-.5564766488918562465935297, -.9148198405846724121600860),
      Complex(-.5564766488918562465935297, +.9148198405846724121600860),
      Complex(-.4345168906815271799687308, -1.045382255856986531461592),
      Complex(-.4345168906815271799687308, +1.045382255856986531461592),
      Complex(-.2640041595834031147954813, -1.192762031948052470183960),
      Complex(-.2640041595834031147954813, +1.192762031948052470183960);
    break;
    
  case 22:
    tempPoles << Complex(-.9058702269930872551848625, -52774908289999045189007100.0e-27),
      Complex(-.9058702269930872551848625, +52774908289999045189007100.0e-27),
      Complex(-.8972983138153530955952835, -.1584351912289865608659759),
      Complex(-.8972983138153530955952835, +.1584351912289865608659759),
      Complex(-.8799661455640176154025352, -.2644363039201535049656450),
      Complex(-.8799661455640176154025352, +.2644363039201535049656450),
      Complex(-.8534754036851687233084587, -.3710389319482319823405321),
      Complex(-.8534754036851687233084587, +.3710389319482319823405321),
      Complex(-.8171682088462720394344996, -.4785619492202780899653575),
      Complex(-.8171682088462720394344996, +.4785619492202780899653575),
      Complex(-.7700332930556816872932937, -.5874255426351153211965601),
      Complex(-.7700332930556816872932937, +.5874255426351153211965601),
      Complex(-.7105305456418785989070935, -.6982266265924524000098548),
      Complex(-.7105305456418785989070935, +.6982266265924524000098548),
      Complex(-.6362427683267827226840153, -.8118875040246347267248508),
      Complex(-.6362427683267827226840153, +.8118875040246347267248508),
      Complex(-.5430983056306302779658129, -.9299947824439872998916657),
      Complex(-.5430983056306302779658129, +.9299947824439872998916657),
      Complex(-.4232528745642628461715044, -1.055755605227545931204656),
      Complex(-.4232528745642628461715044, +1.055755605227545931204656),
      Complex(-.2566376987939318038016012, -1.197982433555213008346532),
      Complex(-.2566376987939318038016012, +1.197982433555213008346532);
    break;
    
  case 23:
    tempPoles << Complex(-.9066732476324988168207439, 0),
      Complex(-.9027564979912504609412993, -.1010534335314045013252480),
      Complex(-.9027564979912504609412993, +.1010534335314045013252480),
      Complex(-.8909283242471251458653994, -.2023024699381223418195228),
      Complex(-.8909283242471251458653994, +.2023024699381223418195228),
      Complex(-.8709469395587416239596874, -.3039581993950041588888925),
      Complex(-.8709469395587416239596874, +.3039581993950041588888925),
      Complex(-.8423805948021127057054288, -.4062657948237602726779246),
      Complex(-.8423805948021127057054288, +.4062657948237602726779246),
      Complex(-.8045561642053176205623187, -.5095305912227258268309528),
      Complex(-.8045561642053176205623187, +.5095305912227258268309528),
      Complex(-.7564660146829880581478138, -.6141594859476032127216463),
      Complex(-.7564660146829880581478138, +.6141594859476032127216463),
      Complex(-.6965966033912705387505040, -.7207341374753046970247055),
      Complex(-.6965966033912705387505040, +.7207341374753046970247055),
      Complex(-.6225903228771341778273152, -.8301558302812980678845563),
      Complex(-.6225903228771341778273152, +.8301558302812980678845563),
      Complex(-.5304922463810191698502226, -.9439760364018300083750242),
      Complex(-.5304922463810191698502226, +.9439760364018300083750242),
      Complex(-.4126986617510148836149955, -1.065328794475513585531053),
      Complex(-.4126986617510148836149955, +1.065328794475513585531053),
      Complex(-.2497697202208956030229911, -1.202813187870697831365338),
      Complex(-.2497697202208956030229911, +1.202813187870697831365338);
    break;
    
  case 24:
    tempPoles << Complex(-.9055312363372773709269407, -48440066540478700874836350.0e-27),
      Complex(-.9055312363372773709269407, +48440066540478700874836350.0e-27),
      Complex(-.8983105104397872954053307, -.1454056133873610120105857),
      Complex(-.8983105104397872954053307, +.1454056133873610120105857),
      Complex(-.8837358034555706623131950, -.2426335234401383076544239),
      Complex(-.8837358034555706623131950, +.2426335234401383076544239),
      Complex(-.8615278304016353651120610, -.3403202112618624773397257),
      Complex(-.8615278304016353651120610, +.3403202112618624773397257),
      Complex(-.8312326466813240652679563, -.4386985933597305434577492),
      Complex(-.8312326466813240652679563, +.4386985933597305434577492),
      Complex(-.7921695462343492518845446, -.5380628490968016700338001),
      Complex(-.7921695462343492518845446, +.5380628490968016700338001),
      Complex(-.7433392285088529449175873, -.6388084216222567930378296),
      Complex(-.7433392285088529449175873, +.6388084216222567930378296),
      Complex(-.6832565803536521302816011, -.7415032695091650806797753),
      Complex(-.6832565803536521302816011, +.7415032695091650806797753),
      Complex(-.6096221567378335562589532, -.8470292433077202380020454),
      Complex(-.6096221567378335562589532, +.8470292433077202380020454),
      Complex(-.5185914574820317343536707, -.9569048385259054576937721),
      Complex(-.5185914574820317343536707, +.9569048385259054576937721),
      Complex(-.4027853855197518014786978, -1.074195196518674765143729),
      Complex(-.4027853855197518014786978, +1.074195196518674765143729),
      Complex(-.2433481337524869675825448, -1.207298683731972524975429),
      Complex(-.2433481337524869675825448, +1.207298683731972524975429);
    break;
    
  case 25:
    tempPoles << Complex(-.9062073871811708652496104, 0),
      Complex(-.9028833390228020537142561, -93077131185102967450643820.0e-27),
      Complex(-.9028833390228020537142561, +93077131185102967450643820.0e-27),
      Complex(-.8928551459883548836774529, -.1863068969804300712287138),
      Complex(-.8928551459883548836774529, +.1863068969804300712287138),
      Complex(-.8759497989677857803656239, -.2798521321771408719327250),
      Complex(-.8759497989677857803656239, +.2798521321771408719327250),
      Complex(-.8518616886554019782346493, -.3738977875907595009446142),
      Complex(-.8518616886554019782346493, +.3738977875907595009446142),
      Complex(-.8201226043936880253962552, -.4686668574656966589020580),
      Complex(-.8201226043936880253962552, +.4686668574656966589020580),
      Complex(-.7800496278186497225905443, -.5644441210349710332887354),
      Complex(-.7800496278186497225905443, +.5644441210349710332887354),
      Complex(-.7306549271849967721596735, -.6616149647357748681460822),
      Complex(-.7306549271849967721596735, +.6616149647357748681460822),
      Complex(-.6704827128029559528610523, -.7607348858167839877987008),
      Complex(-.6704827128029559528610523, +.7607348858167839877987008),
      Complex(-.5972898661335557242320528, -.8626676330388028512598538),
      Complex(-.5972898661335557242320528, +.8626676330388028512598538),
      Complex(-.5073362861078468845461362, -.9689006305344868494672405),
      Complex(-.5073362861078468845461362, +.9689006305344868494672405),
      Complex(-.3934529878191079606023847, -1.082433927173831581956863),
      Complex(-.3934529878191079606023847, +1.082433927173831581956863),
      Complex(-.2373280669322028974199184, -1.211476658382565356579418),
      Complex(-.2373280669322028974199184, +1.211476658382565356579418);
    break;
    
  default:
    // Throw error ValueError the Bessel filter is only supported up to order 25
    break;
  }
  
  (*poles) = tempPoles.rowwise().replicate( channels );
  
}

void coeffsToZpk(const MatrixXR&  b, const MatrixXR&  a, MatrixXC* zeros, MatrixXC* poles, Real* gain){
  (*gain) = b(0, 0);
  MatrixXR bTemp = b;
  bTemp /= b(0, 0);
  roots(bTemp, zeros);
  roots(a, poles);
}

void zpkToCoeffs(const MatrixXC& zeros, const MatrixXC& poles, Real gain, MatrixXC*  b, MatrixXC*  a){
  poly( zeros, b );
  (*b) *= gain;
  
  poly( poles, a );
}

void lowPassToLowPass(const MatrixXC& b, const MatrixXC& a, Real freq, MatrixXC*  bout, MatrixXC*  aout) {
  const int asize = a.cols();
  const int bsize = b.cols();
  
  const int rows = a.rows();
  
  const int maxsize = max(asize, bsize);    
  
  *aout = a;
  *bout = b;
  
  MatrixXR pwo;
  range(maxsize-1, -1, maxsize, rows, &pwo);
  
  pwo = pwo.cwise().expN( freq );
  
  int start1 = max(bsize - asize, 0);
  int start2 = max(asize - bsize, 0);
  
  for ( int i = 0; i < bsize; i++ ) {
    (*bout).col(i).cwise() *= pwo.col( start2 + i ).cwise().inverse() * pwo.col( start1 );
  }
  
  for ( int i = 0; i < asize; i++ ) {
    (*aout).col(i).cwise() *= pwo.col( start1 + i ).cwise().inverse() * pwo.col( start1 );
  }
  
  normalize((*bout), (*aout));
}

void lowPassToHighPass(const MatrixXC& b, const MatrixXC& a, Real freq, MatrixXC*  bout, MatrixXC*  aout) {
  const int asize = a.cols();
  const int bsize = b.cols();
  
  const int rows = a.rows();
  
  const int maxsize = max(asize, bsize);    
  
  (*aout) = MatrixXC::Zero(rows, maxsize);
  (*bout) = MatrixXC::Zero(rows, maxsize);
  
  (*aout).block(0, 0, rows, asize) = a.rowwise().reverse();
  (*bout).block(0, 0, rows, bsize) = b.rowwise().reverse();
  
  MatrixXR pwo;
  range(0, maxsize, maxsize, rows, &pwo);
  
  pwo = pwo.cwise().expN( freq );
  
  (*aout) = (*aout).cwise() * pwo;
  (*bout) = (*bout).cwise() * pwo;
  
  normalize((*bout), (*aout));
}

void lowPassToBandPass(const MatrixXC& b, const MatrixXC& a, Real freq, Real bandwidth, MatrixXC*  bout, MatrixXC*  aout) { 
  const int asize = a.cols();
  const int bsize = b.cols();
  
  const int rows = a.rows();
  
  const int maxsize = max(asize - 1, bsize - 1);
  
  const int maxasize = maxsize + asize - 1;
  const int maxbsize = maxsize + bsize - 1;
  
  const Real freqSq = freq * freq;
  
  (*aout) = MatrixXC::Zero(rows, maxasize + 1);
  (*bout) = MatrixXC::Zero(rows, maxbsize + 1);
  
  for ( int j = 0; j < (maxbsize + 1); j++ ) {
    MatrixXC val = MatrixXC::Zero(rows, 1);
    
    for ( int i = 0; i < bsize; i++ ) {
      for ( int k = 0; k < (i + 1); k++ ) {
        if ( maxsize - i + 2 * k == j ) {
          val += (Real)combination(i, k) * b.col(bsize - 1 - i) * pow(freqSq, (Real)(i-k)) / pow(bandwidth, (Real)i);
        }
      }
    }
    
    (*bout).col(maxbsize - j) = val;
  }
  
  for ( int j = 0; j < (maxasize + 1); j++ ) {
    MatrixXC val = MatrixXC::Zero(rows, 1);
    
    for ( int i = 0; i < asize; i++ ) {
      for ( int k = 0; k < (i + 1); k++ ) {
        if ( maxsize - i + 2 * k == j ) {
          val += (Real)combination(i, k) * a.col(asize - 1 - i) * pow(freqSq, (Real)(i-k)) / pow(bandwidth, (Real)i);
        }
      }
    }
    
    (*aout).col(maxasize - j) = val;
  }
  
  normalize((*bout), (*aout));
}


void lowPassToBandStop(const MatrixXC& b, const MatrixXC& a, Real freq, Real bandwidth, MatrixXC*  bout, MatrixXC*  aout) { 
  const int asize = a.cols();
  const int bsize = b.cols();
  
  const int rows = a.rows();
  
  const int maxsize = max(asize - 1, bsize - 1);
  
  const int maxasize = 2 * maxsize;
  const int maxbsize = 2 * maxsize;
  
  const Real freqSq = freq * freq;
  
  (*aout) = MatrixXC::Zero(rows, maxasize + 1);
  (*bout) = MatrixXC::Zero(rows, maxbsize + 1);
  
  for ( int j = 0; j < (maxbsize + 1); j++ ) {
    MatrixXC val = MatrixXC::Zero(rows, 1);
    
    for ( int i = 0; i < bsize; i++ ) {
      for ( int k = 0; k < (maxsize - i + 1); k++ ) {
        if ( i + 2 * k == j ) {
          val += (Real)combination(maxsize - i, k) * b.col(bsize - 1 - i) * pow(freqSq, (Real)(maxsize - i - k)) * pow(bandwidth, (Real)i);
        }
      }
    }
    
    (*bout).col(maxbsize - j) = val;
  }
  
  for ( int j = 0; j < (maxasize + 1); j++ ) {
    MatrixXC val = MatrixXC::Zero(rows, 1);
    
    for ( int i = 0; i < asize; i++ ) {
      for ( int k = 0; k < (maxsize - i + 1); k++ ) {
        if ( i + 2 * k == j ) {
          val += (Real)combination(maxsize - i, k) * a.col(asize - 1 - i) * pow(freqSq, (Real)(maxsize - i - k)) * pow(bandwidth, (Real)i);
        }
      }
    }
    
    (*aout).col(maxasize - j) = val;
  }
  
  
  normalize((*bout), (*aout));
}


void normalize(MatrixXC& b, MatrixXC& a) {
  
  for (int i = 0; i < b.cols(); i++ ) {
    b.col(i).cwise() /= a.col(0);
  }
  
  for (int i = a.cols()-1; i >= 0; i-- ) {
    a.col(i).cwise() /= a.col(0);
  }
}

void bilinear(const MatrixXC& b, const MatrixXC& a, Real fs, MatrixXR*  bout, MatrixXR*  aout) {
  const int asize = a.cols();
  const int bsize = b.cols();
  const int maxsize = max(asize, bsize);
  
  const int rows = a.rows();
  
  (*aout).resize(rows, maxsize);
  (*bout).resize(rows, maxsize);
  
  MatrixXC val;
  
  for ( int j = 0; j < maxsize; j++ ) {
    val = MatrixXC::Zero(rows, 1);
    
    for ( int i = 0; i < bsize; i++ ) {
      for ( int k = 0; k < i + 1; k++ ) {
        for ( int l = 0; l < (maxsize - i); l++ ) {
          
          if((k + l) == j)
            val += combination(i, k) * combination(maxsize - i - 1, l) * b.col(bsize - 1 - i) * pow(2*fs, (Real)i) * pow((Real)-1, (Real)k);
        }
      }
    }
    
    (*bout).col(j) = val.real();
  }
  
  for ( int j = 0; j < maxsize; j++ ) {
    val = MatrixXC::Zero(rows, 1);
    
    for ( int i = 0; i < asize; i++ ) {
      for ( int k = 0; k < i + 1; k++ ) {
        for ( int l = 0; l < (maxsize - i); l++ ) {
          
          if((k + l) == j)
            val += combination(i, k) * combination(maxsize - i - 1, l) * a.col(asize - 1 - i) * pow((Real)2.0*fs, (Real)i) * pow((Real)-1, (Real)k);
        }
      }
    }
    
    (*aout).col(j) = val.real();
  }
  
}

void bilinear(const MatrixXC& b, const MatrixXC& a, MatrixXR*  bout, MatrixXR*  aout) {
  bilinear(b, a, 1.0, bout, aout);
}
