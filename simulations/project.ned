//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with this program.  If not, see http://www.gnu.org/licenses/.
// 

package teste_mapa.simulations;

import inet.common.geometry.common.OsgGeographicCoordinateSystem;
import inet.environment.common.PhysicalEnvironment;
import inet.networklayer.configurator.ipv4.Ipv4NetworkConfigurator;
import inet.node.inet.AdhocHost;
import inet.physicallayer.ieee80211.packetlevel.Ieee80211ScalarRadioMedium;
import inet.visualizer.integrated.IntegratedMultiVisualizer;
import inet.visualizer.integrated.IntegratedVisualizer;
import inet.visualizer.base.MediumVisualizerBase; //carregar o modulo de visualizacao osg
import inet.visualizer.contract.IMediumVisualizer;
import inet.physicallayer.common.packetlevel.RadioMedium;

network Network
{
    @display("bgb=3574.3176,1851.055");
    @statistic[receptionPower](source="receptionMinSignalPower(radioMedium.signalArrivalStarted)"; record=vector);
    @figure[receptionPower](type=text; pos=5,35; text="potencia: N/A");
    submodules:
        source: AdhocHost {
            @display("p=1475.4226,832.5875");
        }
        destination: AdhocHost {
            @display("p=1475.4226,832.5875;i=misc/sensor");
        }
        configurator: Ipv4NetworkConfigurator {
            @display("p=153.44875,86.45");
        }
        radioMedium: Ieee80211ScalarRadioMedium {
            @display("p=146.965,200.99625");
        }
        physicalEnvironment: PhysicalEnvironment {
            @display("p=146.965,354.44498");
        }
        visualizer: IntegratedVisualizer {
            @display("p=140.48125,523.02246");
        }
        //integratedMultiVisualizer: IntegratedMultiVisualizer {
          //  @display("p=95.095,682.95496");
        //}
        coordinateSystem: OsgGeographicCoordinateSystem {
            parameters:
                @display("p=100,150");
        }
      
}

