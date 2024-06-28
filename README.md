# ACMSI
**The shiny-based software for Fragment Analysis(CE), especially for MSI analysis<br><br>**
The `Windows` version for directly used:<br>
[click here to download ACMSI](https://szfile.haplox.net:7071/f/108c6e555eb7495ea72c/?dl=1 "windows version")<br>
(`Linux or Mac` can run **app.R** with **MSIanalysis.py** to build a server)
<br>
<table border="1.5">
    <tr>
        <th rowspan="3">Software</th>
        <th rowspan="3">Size Calling (322)</th>
        <th colspan="4" align="center">MSI Classification</th>
    </tr>
    <tr>
        <th colspan="2">MSI markers (725)</th>
        <th colspan="2">Samples (126)</th>
    </tr>
    <tr>
        <td>sensitivity</td>
        <td>specificity</td>
        <td>sensitivity</td>
        <td>specificity</td>
    </tr>
    <tr>
        <td><b>ACMSI</b></td>
        <td>99.07%</td>
        <td>96.30%</td>
        <td>99.64%</td>
        <td>100%</td>
        <td>100%</td>
    </tr>
    <tr>
        <td><b>GeneMarker</b></td>
        <td>97.52%</td>
        <td>-</td>
        <td>-</td>
        <td>-</td>
        <td>-</td>
    </tr>
    <tr>
        <td><b>PeakScanner</b></td>
        <td>97.52% (0.93%)</td>
        <td>-</td>
        <td>-</td>
        <td>-</td>
        <td>-</td>
    </tr>
    <tr>
        <td><b>GeneMapper</b></td>
        <td>95.96%</td>
        <td>-</td>
        <td>-</td>
        <td>-</td>
        <td>-</td>
    </tr>
    <tr>
        <td><b>coffalyser.NET</b></td>
        <td>46.27%</td>
        <td>-</td>
        <td>-</td>
        <td>-</td>
        <td>-</td>
    </tr>
</table>

<br>

**The overall workflow shows:<br>**
![The workflow of ACMSI](https://github.com/CrazyJayyy/ACMSI/assets/173884768/61208972-613e-407d-b568-9404fc6b9295)

**The result of MSI classification:<br>**
The `red` portion represents **peaks shared** between tumor and normal tissues<br> 
The `purple` portion represents **elevated peaks** in tumor tissue compared to normal tissue<br> 
The `black` portion represents **new peaks** in tumor tissue<br>
![图片](https://github.com/OpenGene/ACMSI/assets/173884768/40169efe-98d4-496e-b347-be3ad61d27d8)

**Followings show the UI of 3 main modules of ACMSI:<br>**

**<br>1. Plot of Fragment Analysis**
![2024-06-14-14-31-10](https://github.com/OpenGene/ACMSI/assets/173884768/d1b857de-3b7b-49fb-b0d7-a1be49def105)

**<br>2. Calibration of Size Calling**
![2024-06-20-09-41-12](https://github.com/OpenGene/ACMSI/assets/173884768/e75febdd-60e0-41a7-bef6-a4db6a20fc7e)

**<br>3. Automated Classification of MSI**
![2024-06-14-14-40-37](https://github.com/OpenGene/ACMSI/assets/173884768/4a0becc5-48e4-4996-a935-e08016c6399f)

