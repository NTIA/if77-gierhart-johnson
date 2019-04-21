# IF-77 Electromagnetic Wave Propagation Model (Gierhart-Johnson)

In 1973, ITS released the IF-73 (ITS-FAA-1973) air/ground propagation model developed for the Federal Aviation Administration (FAA). IF-73 evolved into the IF-77 (ITS-FAA-1977) model, which is applicable to air/air, air/ground, air/satellite, ground/ground, and ground/satellite paths. The IF-77 was incorporated into a number of FORTRAN programs that are useful in estimating the service coverage of radio systems. These programs may be used to obtain a wide variety of computer-generated microfilm plots.

The IF-77 propagation model is applicable to air/ground, air/air, ground/satellite, and air/satellite paths. It can also be used for ground/ground paths that are line-of-sight or smooth earth. Model applications are restricted to telecommunication systems operating at radio frequencies from about 0.1 to 20 GHz with antenna heights greater than 0.5 m. In addition, radio-horizon elevations must be less than the elevation of the higher antenna. The radio horizon for the higher antenna is taken either as a common horizon with the lower antenna or as a smooth earth horizon with the same elevation as the lower antenna effective reflecting plane.

The work was performed under contract to the U.S. Department of Transportation (DOT), contract number DTFA01-82-Y-10539 and published in the report:

* M.E. Johnson and G.D. Gierhart, [The IF-77 Electromagnetic Wave Propagation Model](https://www.its.bldrdoc.gov/publications/details.aspx?pub=2524), NTIA Sponsor Report FAA-ES-83/3, September 1983.

The IF-77 software is available for download as a "zip" file. The download contains software developed by NTIA. NTIA does not make any warranty of any kind, express, implied or statutory, including, without limitation, the implied warranty of merchantability, fitness for a particular purpose, non-infringement and data accuracy. NTIA does not warrant or make any representations regarding the use of the software or the results thereof, including but not limited to the correctness, accuracy, reliability or usefulness of the software or the results. You can use, copy, modify, and redistribute the NTIA-developed software upon your acceptance of these terms and conditions and upon your express agreement to provide appropriate acknowledgments of NTIA's ownership of and development of the software by keeping this exact text present in any copied or derivative works. The "zip" file is available for download in the [Releases](https://github.com/NTIA/if77-legacy/releases) page of this code repository.

The zip package contains:

* ata.exe: An executable file for ATA that runs on a PC (tested on Windows 7)
* atoa.in: Sample input file for ATA. Run with: ata < atoa.in > atoa.out
* atoa.out: Sample output file for ATA.
* cards.txt: Text description of data input cards.

Extensive comparisons of IF-77 predictions with measured data were made, and an atlas of basic transmission loss predictions was generated using the model. The following publications provide additional information:

* M.E. Johnson and G.D. Gierhart, [Applications Guide for Propagation and Interference Analysis Computer Programs (0.1 to 20 GHz)](https://www.its.bldrdoc.gov/publications/2516.aspx), NTIA Sponsor Report FAA-RD-77-60, March 1978
* M. E. Johnson and G. D. Gierhart, "[Aerospace Propagation Prediction Capabilities Associated with the IF-77 Model](https://www.its.bldrdoc.gov/publications/2683.aspx)," in Conference Proceedings: Operational Modelling of the Aerospace Propagation Environment, North Atlantic Treaty Organization (NATO), Advisory Group for Aerospace Research and Development (AGARD), Neuilly-sur-Seine (France), vol. 2, no. 238, November 1978
* M.E. Johnson and G.D. Gierhart, [Comparison of Measured Data with IF-77 Propagation Model Predictions](https://www.its.bldrdoc.gov/publications/2518.aspx), NTIA Sponsor Report FAA-RD-79-9, August 1979
* M.E. Johnson and G.D. Gierhart, [An Atlas of Basic Transmission Loss For 0.125 to 15.5 GHz](https://www.its.bldrdoc.gov/publications/2520.aspx), NTIA Sponsor Report FAA-RD-80-1, August 1980
* Weiner, M., [Use of the Longley-Rice and Johnson-Gierhart Tropospheric Radio Propagation Programs: 0.02-20 GHz](https://www.its.bldrdoc.gov/media/35805/Weiner.UseL-R&J-G.pdf), IEEE Journal on Selected Areas in Communications, vol.4, no.2, pp. 297-307, March 1986

For assistance with downloads, contact info@its.bldrdoc.gov. For technical questions about IF-77, contact Paul McKenna, (303) 497-3474, pmckenna@ntia.gov.
