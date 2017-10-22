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
#include <QDebug>

PhantomReconstruction::PhantomReconstruction(QWidget *parent) : QWidget(parent),
	m_imageLabel(new QLabel)
{
	setAutoFillBackground(true);

	QPalette palette;
	palette.setColor(QPalette::Background, Qt::black);
	setPalette(palette);

	m_imageLabel->setSizePolicy(QSizePolicy::Ignored, QSizePolicy::Ignored);
	m_imageLabel->setScaledContents(true);

//	QScopedPointer<QVBoxLayout> newLayout(new QVBoxLayout);
//	newLayout->addWidget(m_imageLabel);
//	setLayout(newLayout.data());

//	m_phantom = std::make_shared<Phantom>(fieldOfView.toStdVector());
	std::vector<float> fieldOfView = {28,28};
	m_phantom2D = std::make_shared<Phantom>(fieldOfView);

	QImage image;
	m_imageLabel->setPixmap(QPixmap::fromImage(image));
}

void PhantomReconstruction::reconstruct(Trajectory* trajectory)
{
	if(!isEnabled())
		return;

	std::vector<float> fieldOfView(std::begin(trajectory->fieldOfView), std::end(trajectory->fieldOfView));
//	auto phantom = std::unique_ptr<Phantom>(new Phantom(fieldOfView));

	std::vector<int> acquisitionSize =
	{trajectory->readoutPoints, trajectory->readouts};

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
					kSpaceData.setSignalValue(m, m_phantom2D->fourierDomainSignal(k[0], k[1]));
					break;
				case 3:
//				kSpaceData.setSignalValue(m, m_phantom3D.fourierDomainSignal(k[0], k[1], k[2]));
					break;
				default:
					break;
			}
		}
	}
	Gridding gridding(trajectory);
	MRdata* imageReconstructed = gridding.kSpaceToImage(kSpaceData)[0];
	imageReconstructed->writeToOctave("phantom.txt");
	std::vector<int> imageDimensions = imageReconstructed->dimensions();
	QVector<float> imageMagnitude;
	float maxMagnitude = 0;
	foreach(complexFloat value, imageReconstructed->signal())
	{
		float magnitude = fabs(value);
		maxMagnitude = qMax(maxMagnitude, magnitude);
		imageMagnitude.append(magnitude);
	}

	QVector<unsigned char> imageData(imageMagnitude.size());
	QImage image(imageDimensions[0], imageDimensions[1], QImage::Format_Grayscale8);
	for(int n=0; n<imageMagnitude.size(); n++)
	{
		int y = n%imageDimensions[0];
		int x = n/imageDimensions[0];
		int value = qRound(imageMagnitude[n]/maxMagnitude*255);
		QRgb color = qRgb(value, value, value);
		image.setPixel(x, y, color);
	}
//	QByteArray imageArray(imageData.constData(), imageMagnitude.size());

	//image.save("phantom.jpg");

	m_imageLabel->setPixmap(QPixmap::fromImage(image));
	m_imageLabel->adjustSize();

	update();
}

void PhantomReconstruction::paintEvent(QPaintEvent *event)
{
	Q_UNUSED(event);

	QPainter painter(this);
	QRect viewport = painter.viewport();
	QSize imageSize = m_imageLabel->pixmap()->size();
	imageSize.scale(viewport.size(), Qt::KeepAspectRatio);
	painter.setViewport(viewport.x(), viewport.y(), imageSize.width(), imageSize.height());
	painter.setWindow(m_imageLabel->pixmap()->rect());
	painter.drawPixmap(0, 0, *m_imageLabel->pixmap());
}
