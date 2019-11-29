# IF-77 Air-to-Ground Propagation Model (Gierhart-Johnson)

## Applicability of the Propagation Model

The IF-77 radio wave propagation model is applicable to air/ground, air/air, ground/satellite, and air/satellite paths. It can also be used for ground/ground paths that are line-of-sight or smooth earth. Model applications are restricted to telecommunication systems operating at radio frequencies from about 0.1 to 20 GHz with antenna heights greater than 0.5 m. In addition, radio-horizon elevations must be less than the elevation of the higher antenna. The radio horizon for the higher antenna is taken either as a common horizon with the lower antenna or as a smooth earth horizon with the same elevation as the lower antenna effective reflecting plane.

## History of the Software Implementation

In 1973, ITS released the IF-73 (ITS-FAA-1973) air/ground propagation model developed for the Federal Aviation Administration (FAA). IF-73 evolved into the IF-77 (ITS-FAA-1977) model. The IF-77 was incorporated into a number of FORTRAN programs that are useful in estimating the service coverage of radio systems. These programs were used to obtain a wide variety of computer-generated microfilm plots. Extensive comparisons of IF-77 predictions with measured data were made, and an atlas of basic transmission loss predictions was generated using the model. The work was performed under contract to the U.S. Department of Transportation (DOT), contract number DTFA01-82-Y-10539 and published in 1983.

## Releases

A list of all available software releases can be viewed [here](https://github.com/NTIA/if77-gierhart-johnson/releases). The original IF-77 software written by NTIA is Release v1.0. The release package contains:

* ata.exe: An executable file for ATA that runs on a PC (tested on Windows 7)
* atoa.in: Sample input file for ATA. Run with: `ata.exe < atoa.in > atoa.out`
* atoa.out: Sample output file for ATA.
* cards.txt: Text description of data input cards.

## References

* M.E. Johnson and G.D. Gierhart, [The IF-77 Electromagnetic Wave Propagation Model](https://www.its.bldrdoc.gov/publications/2524.aspx), NTIA Sponsor Report FAA-ES-83/3, September 1983.</a>
* M.E. Johnson and G.D. Gierhart, [Applications Guide for Propagation and Interference Analysis Computer Programs (0.1 to 20 GHz)](https://www.its.bldrdoc.gov/publications/2516.aspx), NTIA Sponsor Report FAA-RD-77-60, March 1978
* M. E. Johnson and G. D. Gierhart, "[Aerospace Propagation Prediction Capabilities Associated with the IF-77 Model](https://www.its.bldrdoc.gov/publications/2683.aspx)," in Conference Proceedings: Operational Modelling of the Aerospace Propagation Environment, North Atlantic Treaty Organization (NATO), Advisory Group for Aerospace Research and Development (AGARD), Neuilly-sur-Seine (France), vol. 2, no. 238, November 1978
* M.E. Johnson and G.D. Gierhart, [Comparison of Measured Data with IF-77 Propagation Model Predictions](https://www.its.bldrdoc.gov/publications/2518.aspx), NTIA Sponsor Report FAA-RD-79-9, August 1979
* M.E. Johnson and G.D. Gierhart, [An Atlas of Basic Transmission Loss For 0.125 to 15.5 GHz](https://www.its.bldrdoc.gov/publications/2520.aspx), NTIA Sponsor Report FAA-RD-80-1, August 1980
* Weiner, M., [Use of the Longley-Rice and Johnson-Gierhart Tropospheric Radio Propagation Programs: 0.02-20 GHz](https://www.its.bldrdoc.gov/media/35805/Weiner.UseL-R&J-G.pdf), IEEE Journal on Selected Areas in Communications, vol.4, no.2, pp. 297-307, March 1986

## Contact

For technical questions about IF-77, contact Paul McKenna, (303) 497-3474, pmckenna@ntia.gov.
