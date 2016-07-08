package cpm_boundary;

import java.awt.*;

import javax.swing.*;

import org.jfree.chart.*;
import org.jfree.chart.plot.*;
import org.jfree.data.xy.*;

@SuppressWarnings("serial")
public class EnergyPanel extends JPanel implements DataListener {
	
	private CellPottsModel model = null;
	private XYSeries series = null;
	private XYSeriesCollection data = null;
	private JFreeChart chart = null;
	private ChartPanel chartPanel = null;
	
	public void setModel(CellPottsModel model){
		if (this.model != null){
			this.model.removeDataListener(this);
		}
		this.model = model;
		this.model.addDataListener(this);
		init();
	}
	
	public void init(){
		if (chartPanel != null){
			this.remove(chartPanel);
		}
		series = new XYSeries("Energy data");
		data = new XYSeriesCollection();
		data.addSeries(series);
		chart = ChartFactory.createXYLineChart(
				"Energy data",
				"Time (MCS)",
				"Energy (kT)",
				data,
				PlotOrientation.VERTICAL,
				true,
				true,
				false);
		chartPanel = new ChartPanel(chart);
		this.setLayout(new BorderLayout());
		this.add(chartPanel, BorderLayout.CENTER);
		this.validate();
	}

	@Override
	public void update(CellPottsModel model, int time) {
		series.add(time, model.calculateRoughness());
	}
}
