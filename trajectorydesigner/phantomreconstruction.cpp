#include "phantomreconstruction.h"

#include "phantom.h"
#include "mrdata.h"
#include "gridding.h"

extern "C"
{
#include "trajectory.h"
}

#include <QVBoxLayout>
#include <QLabel>
#include <QImage>
#include <QPainter>

PhantomReconstruction::PhantomReconstruction(QWidget *parent) : QWidget(parent),
	m_imageLabel(new QLabel)
{
	setAutoFillBackground(true);

	QPalette palette;
	palette.setColor(QPalette::Background, Qt::black);
	setPalette(palette);

	QScopedPointer<QVBoxLayout> newLayout(new QVBoxLayout);
	newLayout->addWidget(m_imageLabel);
	setLayout(newLayout.data());

//	m_phantom = std::make_shared<Phantom>(fieldOfView.toStdVector());
	QImage image;
	m_imageLabel->setPixmap(QPixmap::fromImage(image));
}

void PhantomReconstruction::reconstruct(Trajectory* trajectory)
{
	std::vector<float> fieldOfView(std::begin(trajectory->fieldOfView), std::end(trajectory->fieldOfView));
//	auto phantom = std::unique_ptr<Phantom>(new Phantom(fieldOfView));

	std::vector<int> acquisitionSize =
	{trajectory->readoutPoints, trajectory->readouts};

	Phantom phantom(fieldOfView);
	MRdata kSpaceData(acquisitionSize, trajectory->dimensions);

	for(int n=0; n<trajectory->readoutPoints; n++)
	{
		for(int r=0; r<trajectory->readouts; r++)
		{
			float k[3];
			float densityCompensation;
			trajectoryCoordinates(n, r, trajectory, k, &densityCompensation);
			int m = trajectory->readoutPoints*r + n;
			switch (trajectory->dimensions) {
				case 2:
					kSpaceData.setSignalValue(m, phantom.fourierDomainSignal(k[0], k[1]));
					break;
				case 3:
				kSpaceData.setSignalValue(m, phantom.fourierDomainSignal(k[0], k[1], k[2]));
					break;
				default:
					break;
			}
		}
	}
	Gridding gridding(trajectory);
	MRdata* imageReconstructed = gridding.kSpaceToImage(kSpaceData)[0];
	std::vector<int> imageDimensions = imageReconstructed->dimensions();
	QImage image(imageDimensions[0], imageDimensions[1], QImage::Format_Grayscale8);
	m_imageLabel->setPixmap(QPixmap::fromImage(image));
	m_imageLabel->adjustSize();
}

void PhantomReconstruction::paintEvent(QPaintEvent *event)
{
	Q_UNUSED(event);

	QPainter painter(this);
	painter.drawPixmap(0, 0, *m_imageLabel->pixmap());
}
