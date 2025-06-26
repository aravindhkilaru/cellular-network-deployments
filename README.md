# cellular-network-deployments

**To implement a Cellular Network Layout for various coverage scenarios:**
  ● Urban Micro (UMi).
  ● Urban Macro (UMa).
  ● Sub-Urban Macro (SMa)
  ● Rural Macro (RMa).
  ● Indoor Factory (InF).

**Evaluate the characteristics and performance of Cellular Network Layouts:**
  ● Calculation of Signal-to-Interference-plus-Noise Ratio (SINR),
  ● Calculation of Reference Signal Received Power (RSRP),
  ● Path Loss Analysis,
  ● Connection Density,
  ● Area Traffic Capacity,
  ● Control plane latency calculation,
  ● User plane latency calculation,
  ● Data Rate,
    o User Experienced Data Rate (throughput).
  ● Air Interface Latency,
  ● Average Spectral Efficiency,
  ● 5th Percentile User Spectral Efficiency (CDF): The cumulative distribution function
  representing the spectral efficiency experienced by the users.

You can use the code in `network_layout.py` to generate cellular network layouts honeycomb grids for Urban Micro (UMi), Urban Macro (UMa), Sub-Urban Macro (SMa), and Rural Macro (RMa).
![honeycomb_grid](https://github.com/user-attachments/assets/2b1053c6-3d58-4762-b50c-1596e2fd8e3d)

You can use the code `5g_4g_network_simulation.py` to generate honeycomb grid layouts of cellular networks for various environments, including Urban Micro (UMi), Urban Macro (UMa), Sub-Urban Macro (SMa), and Rural Macro (RMa). This code allows you to evaluate the characteristics and performance of these cellular network layouts by calculating the following metrics: Signal-to-Interference-plus-Noise Ratio (SINR), Reference Signal Received Power (RSRP), path loss, connection density, area traffic capacity, control plane latency, user plane latency, data rate, user experienced data rate (throughput), air interface latency, average spectral efficiency, and the 5th percentile user spectral efficiency, which represents the cumulative distribution function of spectral efficiency experienced by users.

![4G_layout_with_UEs](https://github.com/user-attachments/assets/30ffcdd8-85b3-4379-9848-d5d64ec2d0ec)
![SINR_vs_Distance](https://github.com/user-attachments/assets/c912f892-87b9-4c32-b572-f0bf26bfd7f3)
![RSRP_VS_Distance](https://github.com/user-attachments/assets/5899c2d1-2e2b-4c49-a084-9af953d7ded6)
![PL_VS_Distance](https://github.com/user-attachments/assets/f4c5c72e-c10f-403d-9177-096b3dccf598)
![spectral_efficiency_CDF](https://github.com/user-attachments/assets/c8568e70-a6d3-4dcb-a62e-0f936f704b50)
![CDF_SINR](https://github.com/user-attachments/assets/ad1dfd11-a1ed-43cb-b803-04e029cce1f8)

and
Connection Density: 69.54 UEs/km²
Average Spectral Efficiency: 0.52 bps/Hz
Area Capacity: 36.22 bps/Hz/km²

**_Note: These codes were designed to teach undergraduate students about cellular network layouts._**
