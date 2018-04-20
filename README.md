# kmplot
Plot the Kaplan-Meier estimation of the survival function</br>
Survival times are data that measure follow-up time from a defined
starting point to the occurrence of a given event, for example the time
from the beginning to the end of a remission period or the time from the
diagnosis of a disease to death. Standard statistical techniques cannot
usually be applied because the underlying distribution is rarely Normal
and the data are often "censored". A survival time is described as
censored when there is a follow-up time but the event has not yet
occurred or is not known to have occurred. For example, if remission time
is being studied and the patient is still in remission at the end of the
study, then that patient's remission time would be censored. If a patient
for some reason drops out of a study before the end of the study period,
then that patient's follow-up time would also be considered to be
censored. The survival function S(t) is defined as the probability of
surviving at least to time t. The graph of S(t) against t is called the
survival curve. The Kaplan-Meier method can be used to estimate this
curve from the observed survival times without the assumption of an
underlying probability distribution.

          Created by Giuseppe Cardillo
          giuseppe.cardillo-edta@poste.it

To cite this file, this would be an appropriate format:Curve
Cardillo G. (2008). KMPLOT: Kaplan-Meier estimation of the survival
function.
http://www.mathworks.com/matlabcentral/fileexchange/22293
